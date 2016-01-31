#!/usr/bin/env python2
import os
import sys
import cPickle

from numpy		import *
from numpy.random import *
from tssb		 import *
from alleles	 import *
from util		 import *

import numpy.random

from util2 import *
from params import *
from printo import *

import argparse
import signal
import tempfile
import threading
import traceback
import time
from datetime import datetime

# num_samples: number of MCMC samples
# mh_itr: number of metropolis-hasting iterations
# rand_seed: random seed (initialization). Set to None to choose random seed automatically.
def start_new_run(state_manager, backup_manager, safe_to_exit, run_succeeded, config, ssm_file, cnv_file, top_k_trees_file, clonal_freqs_file, burnin_samples, num_samples, mh_itr, mh_std, write_state_every, write_backups_every, rand_seed, tmp_dir):
	state = {}

	with open('random_seed.txt', 'w') as seedf:
		seedf.write('%s\n' % rand_seed)
	try:
		rand_seed = int(rand_seed)
		state['rand_seed'] = rand_seed
		seed(state['rand_seed'])
	except TypeError:
		# If rand_seed is not provided as command-line arg, it will be None,
		# meaning it will hit this code path. Explicitly avoid calling seed(None)
		# -- though this is currently the equivalent of calling seed() in that it
		# seeds the PRNG with /dev/urandom, the semantics of seed(None) might
		# change in later NumPy versions to always seed to the same state.
		state['rand_seed'] = rand_seed
		seed()

	state['ssm_file'] = ssm_file
	state['cnv_file'] = cnv_file
	state['tmp_dir'] = tmp_dir
	state['top_k_trees_file'] = top_k_trees_file
	state['clonal_freqs_file'] = clonal_freqs_file
	state['write_state_every'] = write_state_every
	state['write_backups_every'] = write_backups_every

	codes, n_ssms, n_cnvs = load_data(state['ssm_file'], state['cnv_file'])
	if len(codes) == 0:
		logmsg('No SSMs or CNVs provided. Exiting.', sys.stderr)
		return
	NTPS = len(codes[0].a) # number of samples / time point
	state['glist'] = [datum.name for datum in codes if len(datum.name)>0]

	# MCMC settings
	state['burnin'] = burnin_samples
	state['num_samples'] = num_samples
	state['dp_alpha'] = 25.0
	state['dp_gamma'] = 1.0
	state['alpha_decay'] = 0.25
	state['top_k'] = 5

	# Metropolis-Hastings settings
	state['mh_burnin'] = 0
	state['mh_itr'] = mh_itr # No. of iterations in metropolis-hastings
	state['mh_std'] = mh_std

	state['cd_llh_traces'] = zeros((state['num_samples'], 1))
	state['burnin_cd_llh_traces'] = zeros((state['burnin'], 1))
	state['working_directory'] = os.getcwd()

	root = alleles(conc=0.1, ntps=NTPS)
	state['tssb'] = TSSB(dp_alpha=state['dp_alpha'], dp_gamma=state['dp_gamma'], alpha_decay=state['alpha_decay'], root_node=root, data=codes)
	# hack...
	if 1:
		depth=0
		state['tssb'].root['sticks'] = vstack([ state['tssb'].root['sticks'], boundbeta(1, state['tssb'].dp_gamma) if depth!=0 else .999])
		state['tssb'].root['children'].append({ 'node': state['tssb'].root['node'].spawn(),
					'main':boundbeta(1.0, (state['tssb'].alpha_decay**(depth+1))*state['tssb'].dp_alpha) if state['tssb'].min_depth <= (depth+1) else 0.0, 
					'sticks' : empty((0,1)),	
					'children' : [] })
		new_node = state['tssb'].root['children'][0]['node']
		for n in range(state['tssb'].num_data):
			state['tssb'].assignments[n].remove_datum(n)
			new_node.add_datum(n)
			state['tssb'].assignments[n] = new_node
	
	for datum in codes:
		datum.tssb = state['tssb']
	
	tree_writer = TreeWriter()
	state_manager.write_initial_state(state)
	logmsg("Starting MCMC run...")
	state['last_iteration'] = -state['burnin'] - 1

	# This will overwrite file if it already exists, which is the desired
	# behaviour for a fresh run.
	with open('mcmc_samples.txt', 'w') as mcmcf:
		mcmcf.write('Iteration\tLLH\tTime\n')

	do_mcmc(state_manager, backup_manager, safe_to_exit, run_succeeded, config, state, tree_writer, codes, n_ssms, n_cnvs, NTPS, tmp_dir)

def resume_existing_run(state_manager, backup_manager, safe_to_exit, run_succeeded, config):
	# If error occurs, restore the backups and try again. Never try more than two
	# times, however -- if the primary file and the backup file both fail, the
	# error is unrecoverable.
	try:
		state = state_manager.load_state()
		tree_writer = TreeWriter(resume_run = True)
	except:
		logmsg('Restoring state failed:', sys.stderr)
		traceback.print_exc()
		logmsg('Restoring from backup and trying again.', sys.stderr)
		backup_manager.restore_backup()

		state = state_manager.load_state()
		tree_writer = TreeWriter(resume_run = True)

	set_state(state['rand_state']) # Restore NumPy's RNG state.
	os.chdir(state['working_directory'])
	codes, n_ssms, n_cnvs = load_data(state['ssm_file'], state['cnv_file'])
	NTPS = len(codes[0].a) # number of samples / time point

	do_mcmc(state_manager, backup_manager, safe_to_exit, run_succeeded, config, state, tree_writer, codes, n_ssms, n_cnvs, NTPS, state['tmp_dir'])

def do_mcmc(state_manager, backup_manager, safe_to_exit, run_succeeded, config, state, tree_writer, codes, n_ssms, n_cnvs, NTPS, tmp_dir_parent):
	start_iter = state['last_iteration'] + 1
	unwritten_trees = []
	mcmc_sample_times = []
	last_mcmc_sample_time = time.time()

	# If --tmp-dir is not specified on the command line, it will by default be
	# None, which will cause mkdtemp() to place this directory under the system's
	# temporary directory. This is the desired behaviour.
	config['tmp_dir'] = tempfile.mkdtemp(prefix='pwgsdataexchange.', dir=tmp_dir_parent)

	for iteration in range(start_iter, state['num_samples']):
		safe_to_exit.set()
		if iteration < 0:
			logmsg(iteration)

		# Referring to tssb as local variable instead of dictionary element is much
		# faster.
		tssb = state['tssb']
		tssb.resample_assignments()
		tssb.cull_tree()
		
		# assign node ids
		wts, nodes = tssb.get_mixture()
		for i, node in enumerate(nodes):
			node.id = i
		
		##################################################
		## some useful info about the tree,
		## used by CNV related computations,
		## to be called only after resampling assignments
		set_node_height(tssb)
		set_path_from_root_to_node(tssb)
		map_datum_to_node(tssb)
		##################################################

		state['mh_acc'] = metropolis(
			tssb,
			state['mh_itr'],
			state['mh_std'],
			state['mh_burnin'],
			n_ssms,
			n_cnvs,
			state['ssm_file'],
			state['cnv_file'],
			state['rand_seed'],
			NTPS,
			config['tmp_dir']
		)
		if float(state['mh_acc']) < 0.08 and state['mh_std'] < 10000:
			state['mh_std'] = state['mh_std']*2.0
			logmsg("Shrinking MH proposals. Now %f" % state['mh_std'])
		if float(state['mh_acc']) > 0.5 and float(state['mh_acc']) < 0.99:
			state['mh_std'] = state['mh_std']/2.0
			logmsg("Growing MH proposals. Now %f" % state['mh_std'])
	
		tssb.resample_sticks()
		tssb.resample_stick_orders()
		tssb.resample_hypers(dp_alpha=True, alpha_decay=True, dp_gamma=True)
 
		last_llh = tssb.complete_data_log_likelihood()
		if iteration >= 0:
			state['cd_llh_traces'][iteration] = last_llh
			if True or mod(iteration, 10) == 0:
				weights, nodes = tssb.get_mixture()
				logmsg(' '.join([str(v) for v in (iteration, len(nodes), state['cd_llh_traces'][iteration], state['mh_acc'], tssb.dp_alpha, tssb.dp_gamma, tssb.alpha_decay)]))
			if argmax(state['cd_llh_traces'][:iteration+1]) == iteration:
				logmsg("%f is best per-data complete data likelihood so far." % (state['cd_llh_traces'][iteration]))
		else:
			state['burnin_cd_llh_traces'][iteration + state['burnin']] = last_llh

		unwritten_trees.append((tssb, iteration, last_llh))
		state['tssb'] = tssb
		state['rand_state'] = get_state()
		state['last_iteration'] = iteration

		new_mcmc_sample_time = time.time()
		mcmc_sample_times.append(new_mcmc_sample_time - last_mcmc_sample_time)
		last_mcmc_sample_time = new_mcmc_sample_time

		# It's not safe to exit while performing file IO, as we don't want
		# trees.zip or the computation state file to become corrupted from an
		# interrupted write.
		safe_to_exit.clear()
		should_write_backup = iteration % state['write_backups_every'] == 0 and iteration != start_iter
		should_write_state = iteration % state['write_state_every'] == 0
		is_last_iteration = (iteration == state['num_samples'] - 1)

		# If backup is scheduled to be written, write both it and full program
		# state regardless of whether we're scheduled to write state this
		# iteration.
		if should_write_backup or should_write_state or is_last_iteration:
			with open('mcmc_samples.txt', 'a') as mcmcf:
				llhs_and_times = [(itr, llh, itr_time) for (tssb, itr, llh), itr_time in zip(unwritten_trees, mcmc_sample_times)]
				llhs_and_times = '\n'.join(['%s\t%s\t%s' % (itr, llh, itr_time) for itr, llh, itr_time in llhs_and_times])
				mcmcf.write(llhs_and_times + '\n')
			tree_writer.write_trees(unwritten_trees)
			state_manager.write_state(state)
			unwritten_trees = []
			mcmc_sample_times = []
			if should_write_backup:
				backup_manager.save_backup()

	safe_to_exit.clear()
	#save the best tree
	print_top_trees(TreeWriter.default_archive_fn, state['top_k_trees_file'], state['top_k'])

	#save clonal frequencies
	freq = dict([(g,[] )for g in state['glist']])
	glist = array(freq.keys(),str)
	glist.shape=(1,len(glist))
	savetxt(state['clonal_freqs_file'] ,vstack((glist, array([freq[g] for g in freq.keys()]).T)), fmt='%s', delimiter=', ')
	state_manager.delete_state_file()

	safe_to_exit.set()
	run_succeeded.set()

def test():
	tssb=cPickle.load(open('ptree'))
	wts,nodes=tssb.get_mixture()	
	for dat in tssb.data:
		print [dat.id, dat.__log_likelihood__(0.5)]

def parse_args():
	parser = argparse.ArgumentParser(
		description='Run PhyloWGS to infer subclonal composition from SSMs and CNVs',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument('-b', '--write-backups-every', dest='write_backups_every', default=100, type=int,
		help='Number of iterations to go between writing backups of program state')
	parser.add_argument('-S', '--write-state-every', dest='write_state_every', default=10, type=int,
		help='Number of iterations between writing program state to disk. Higher values reduce IO burden at the cost of losing progress made if program is interrupted.')
	parser.add_argument('-k', '--top-k-trees', dest='top_k_trees', default='top_k_trees',
		help='Output file to save top-k trees in text format')
	parser.add_argument('-f', '--clonal-freqs', dest='clonal_freqs', default='clonalFrequencies',
		help='Output file to save clonal frequencies')
	parser.add_argument('-B', '--burnin-samples', dest='burnin_samples', default=1000, type=int,
		help='Number of burnin samples')
	parser.add_argument('-s', '--mcmc-samples', dest='mcmc_samples', default=2500, type=int,
		help='Number of MCMC samples')
	parser.add_argument('-i', '--mh-iterations', dest='mh_iterations', default=5000, type=int,
		help='Number of Metropolis-Hastings iterations')
	parser.add_argument('-r', '--random-seed', dest='random_seed', type=int,
		help='Random seed for initializing MCMC sampler')
	parser.add_argument('-t', '--tmp-dir', dest='tmp_dir',
		help='Path to directory for temporary files')
	parser.add_argument('ssm_file',
		help='File listing SSMs (simple somatic mutations, i.e., single nucleotide variants. For proper format, see README.md.')
	parser.add_argument('cnv_file',
		help='File listing CNVs (copy number variations). For proper format, see README.md.')
	args = parser.parse_args()
	return args

def run(safe_to_exit, run_succeeded, config):
	state_manager = StateManager()
	backup_manager = BackupManager([StateManager.default_last_state_fn, TreeWriter.default_archive_fn])

	if state_manager.state_exists():
		logmsg('Resuming existing run. Ignoring command-line parameters.')
		resume_existing_run(state_manager, backup_manager, safe_to_exit, run_succeeded, config)
	else:
		args = parse_args()
		# Ensure input files exist and can be read.
		try:
			ssm_file = open(args.ssm_file)
			cnv_file = open(args.cnv_file)
			ssm_file.close()
			cnv_file.close()
		except IOError as e:
			sys.stderr.write(str(e) + '\n')
			sys.exit(1)

		start_new_run(
			state_manager,
			backup_manager,
			safe_to_exit,
			run_succeeded,
			config,
			args.ssm_file,
			args.cnv_file,
			top_k_trees_file=args.top_k_trees,
			clonal_freqs_file=args.clonal_freqs,
			burnin_samples=args.burnin_samples,
			num_samples=args.mcmc_samples,
			mh_itr=args.mh_iterations,
			mh_std=100,
			write_state_every=args.write_state_every,
			write_backups_every=args.write_backups_every,
			rand_seed=args.random_seed,
			tmp_dir=args.tmp_dir
		)

def remove_tmp_files(tmp_dir):
	if tmp_dir is None:
		return
	tmp_filenames = get_c_fnames(tmp_dir)
	for tmpfn in tmp_filenames:
		try:
			os.remove(tmpfn)
		except OSError:
			pass
	try:
		os.rmdir(tmp_dir)
	except OSError:
		pass

def main():
	# Introducing threading is necessary to allow write operations to complete
	# when interrupts are received. As the interrupt halts execution of the main
	# thread and immediately jumps to the interrupt handler, we must run the
	# PhyloWGS code in a different thread, which clears the safe_to_exit flag
	# when in the midst of a write operation. This way, the main thread is left
	# only to handle the signal, allowing the derived thread to finish its
	# current write operation.
	safe_to_exit = threading.Event()

	# This will allow us to detect whether the run thread exited cleanly or not.
	# This means we can properly report a non-zero exit code if something failed
	# (e.g., something threw an exception). This is necessary because exceptions
	# in the run thread will terminate it, but can't be detected from the main
	# thread. A more robust strategy is here: http://stackoverflow.com/a/2830127.
	# Our strategy should be sufficient for the moment, though.
	run_succeeded = threading.Event()

	# We must know where temporary files are stored from within main() so that we
	# can remove them when we exit. However, we don't know this location until
	# the run thread starts, as when PWGS resumes an existing run, the parent
	# directory for the temporary files is stored in the state pickle file. Thus,
	# the run thread will set this value once it is established.
	#
	# So long as this dictionary is used only as a key-value store for primitive
	# objects, it's thread safe and doesn't require the use of a mutex. See
	# http://effbot.org/pyfaq/what-kinds-of-global-value-mutation-are-thread-safe.htm.
	# If more complex values are stored here, we must introduce a mutex.
	config = {
		'tmp_dir': None
	}

	def sigterm_handler(_signo, _stack_frame):
		logmsg('Signal %s received.' % _signo, sys.stderr)
		safe_to_exit.wait()
		remove_tmp_files(config['tmp_dir'])
		logmsg('Exiting now.')
		# Exit with non-zero to indicate run didn't finish.
		sys.exit(3)

	# SciNet will supposedly send SIGTERM 30 s before hard-killing the process.
	# This gives us time to clean up.
	signal.signal(signal.SIGTERM, sigterm_handler)
	# SIGINT is sent on CTRL-C. We don't want the user to interrupt a write
	# operation by hitting CTRL-C, thereby potentially resulting in corrupted
	# data being written. Permit these operations to finish before exiting.
	signal.signal(signal.SIGINT, sigterm_handler)

	run_thread = threading.Thread(target=run, args=(safe_to_exit, run_succeeded, config))
	# Thread must be a daemon thread, or sys.exit() will wait until the thread
	# finishes execution completely.
	run_thread.daemon = True
	run_thread.start()

	while True:
		if not run_thread.is_alive():
			break
		# I don't fully understand this. At least on the imacvm machine, calling
		# join with no timeout argument doesn't work, as the signal handler does
		# not seem to run until run_thread exits. If, however, I specify *any*
		# timeout, no matter how short or long, the signal handler will run
		# *immediately* when the signal is sent -- i.e., even before the timeout
		# has expired.
		run_thread.join(10)

	remove_tmp_files(config['tmp_dir'])
	if run_succeeded.is_set():
		logmsg('Run succeeded.')
		sys.exit(0)
	else:
		logmsg('Run failed.')
		sys.exit(1)

def logmsg(msg, fd=sys.stdout):
	  print >> fd, '[%s] %s' % (datetime.now(), msg)

if __name__ == "__main__":
	main()
