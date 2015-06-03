#!/usr/bin/env python2
import os
import sys
import time
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
import threading
import signal
import traceback

# num_samples: number of MCMC samples
# mh_itr: number of metropolis-hasting iterations
# rand_seed: random seed (initialization). Set to None to choose random seed automatically.
def start_new_run(state_manager, backup_manager, safe_to_exit, ssm_file, cnv_file, top_k_trees_file, clonal_freqs_file, num_samples, mh_itr, mh_std, write_backups_every, rand_seed):
	state = {}
	state['rand_seed'] = rand_seed
	seed(state['rand_seed'])

	state['ssm_file'] = ssm_file
	state['cnv_file'] = cnv_file
	state['top_k_trees_file'] = top_k_trees_file
	state['clonal_freqs_file'] = clonal_freqs_file
	state['write_backups_every'] = write_backups_every

	codes, n_ssms, n_cnvs = load_data(state['ssm_file'], state['cnv_file'])
	NTPS = len(codes[0].a) # number of samples / time point
	state['glist'] = [datum.name for datum in codes if len(datum.name)>0]

	# MCMC settings
	state['burnin'] = 1000
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
	print("Starting MCMC run...")
	state['last_iteration'] = -state['burnin'] - 1

	do_mcmc(state_manager, backup_manager, safe_to_exit, state, tree_writer, codes, n_ssms, n_cnvs, NTPS)

def resume_existing_run(state_manager, backup_manager, safe_to_exit):
	# If error occurs, restore the backups and try again. Never try more than two
	# times, however -- if the primary file and the backup file both fail, the
	# error is unrecoverable.
	try:
		state = state_manager.load_state()
		tree_writer = TreeWriter(resume_run = True)
	except:
		print('Restoring state failed:')
		traceback.print_exc()
		print('Restoring from backup and trying again.')
		backup_manager.restore_backup()

		state = state_manager.load_state()
		tree_writer = TreeWriter(resume_run = True)

	set_state(state['rand_state']) # Restore NumPy's RNG state.
	os.chdir(state['working_directory'])
	codes, n_ssms, n_cnvs = load_data(state['ssm_file'], state['cnv_file'])
	NTPS = len(codes[0].a) # number of samples / time point

	do_mcmc(state_manager, backup_manager, safe_to_exit, state, tree_writer, codes, n_ssms, n_cnvs, NTPS)

def do_mcmc(state_manager, backup_manager, safe_to_exit, state, tree_writer, codes, n_ssms, n_cnvs, NTPS):
	start_iter = state['last_iteration'] + 1

	for iteration in range(start_iter, state['num_samples']):
		safe_to_exit.set()
		if iteration < 0:
			print iteration

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
		)
		if float(state['mh_acc']) < 0.08 and state['mh_std'] < 10000:
			state['mh_std'] = state['mh_std']*2.0
			print "Shrinking MH proposals. Now %f" % state['mh_std']
		if float(state['mh_acc']) > 0.5 and float(state['mh_acc']) < 0.99:
			state['mh_std'] = state['mh_std']/2.0
			print "Growing MH proposals. Now %f" % state['mh_std']
	
		tssb.resample_sticks()
		tssb.resample_stick_orders()
		tssb.resample_hypers(dp_alpha=True, alpha_decay=True, dp_gamma=True)
 
		if iteration >= 0:
			state['cd_llh_traces'][iteration] = tssb.complete_data_log_likelihood()
			if True or mod(iteration, 10) == 0:
				weights, nodes = tssb.get_mixture()
				print iteration, len(nodes), state['cd_llh_traces'][iteration], state['mh_acc'], tssb.dp_alpha, tssb.dp_gamma, tssb.alpha_decay
			if argmax(state['cd_llh_traces'][:iteration+1]) == iteration:
				print "\t%f is best per-data complete data likelihood so far." % (state['cd_llh_traces'][iteration])

		# It's not safe to exit while performing file IO, as we don't want
		# trees.zip or the computation state file to become corrupted from an
		# interrupted write.
		safe_to_exit.clear()
		if iteration >= 0:
			tree_writer.write_tree(tssb, state['cd_llh_traces'][iteration][0], iteration)
		else:
			tree_writer.write_burnin_tree(tssb, iteration)

		state['tssb'] = tssb
		state['rand_state'] = get_state()
		state['last_iteration'] = iteration
		state_manager.write_state(state)

		if iteration % state['write_backups_every'] == 0 and iteration != start_iter:
			backup_manager.save_backup()

	safe_to_exit.clear()
	#save the best tree
	print_top_trees(tree_writer.archive_fn, state['top_k_trees_file'], state['top_k'])

	#save clonal frequencies
	freq = dict([(g,[] )for g in state['glist']])
	glist = array(freq.keys(),str)
	glist.shape=(1,len(glist))
	savetxt(state['clonal_freqs_file'] ,vstack((glist, array([freq[g] for g in freq.keys()]).T)), fmt='%s', delimiter=', ')
	state_manager.delete_state_file()
	safe_to_exit.set()

def test():
	tssb=cPickle.load(open('ptree'))
	wts,nodes=tssb.get_mixture()	
	for dat in tssb.data:
		print [dat.id, dat.__log_likelihood__(0.5)]

def run(safe_to_exit):
	parser = argparse.ArgumentParser(
		description='Run PhyloWGS to infer subclonal composition from SSMs and CNVs',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument('-b', '--write-backups-every', dest='write_backups_every', default=100, type=int,
		help='Number of iterations to go between writing backups of program state')
	parser.add_argument('-k', '--top-k-trees', dest='top_k_trees', default='top_k_trees',
		help='Output file to save top-k trees in text format')
	parser.add_argument('-f', '--clonal-freqs', dest='clonal_freqs', default='clonalFrequencies',
		help='Output file to save clonal frequencies')
	parser.add_argument('-s', '--mcmc-samples', dest='mcmc_samples', default=2500, type=int,
		help='Number of MCMC samples')
	parser.add_argument('-i', '--mh-iterations', dest='mh_iterations', default=5000, type=int,
		help='Number of Metropolis-Hastings iterations')
	parser.add_argument('-r', '--random-seed', dest='random_seed', default=1, type=int,
		help='Random seed for initializing MCMC sampler')
	parser.add_argument('ssm_file',
		help='File listing SSMs (simple somatic mutations, i.e., single nucleotide variants. For proper format, see README.md.')
	parser.add_argument('cnv_file',
		help='File listing CNVs (copy number variations). For proper format, see README.md.')
	args = parser.parse_args()

	# Ensure input files exist and can be read.
	try:
		ssm_file = open(args.ssm_file)
		cnv_file = open(args.cnv_file)
		ssm_file.close()
		cnv_file.close()
	except IOError as e:
		print(e)
		sys.exit(1)

	state_manager = StateManager()
	backup_manager = BackupManager([StateManager.default_last_state_fn, TreeWriter.default_archive_fn])

	if state_manager.state_exists():
		resume_existing_run(state_manager, backup_manager, safe_to_exit)
	else:
		start_new_run(
			state_manager,
			backup_manager,
			safe_to_exit,
			args.ssm_file,
			args.cnv_file,
			top_k_trees_file=args.top_k_trees,
			clonal_freqs_file=args.clonal_freqs,
			num_samples=args.mcmc_samples,
			mh_itr=args.mh_iterations,
			mh_std=100,
			write_backups_every=args.write_backups_every,
			rand_seed=args.random_seed
		)

def main():
	# Introducing threading is necessary to allow write operations to complete
	# when interrupts are received. As the interrupt halts execution of the main
	# thread and immediately jumps to the interrupt handler, we must run the
	# PhyloWGS code in a different thread, which clears the safe_to_exit flag
	# when in the midst of a write operation. This way, the main thread is left
	# only to handle the signal, allowing the derived thread to finish its
	# current write operation.
	safe_to_exit = threading.Event()

	def sigterm_handler(_signo, _stack_frame):
		safe_to_exit.wait()
		# Exit with non-zero to indicate run didn't finish.
		sys.exit(3)
	# SciNet will supposedly send SIGTERM 30 s before hard-killing the process.
	# This gives us time to clean up.
	signal.signal(signal.SIGTERM, sigterm_handler)
	# SIGINT is sent on CTRL-C. We don't want the user to interrupt a write
	# operation by hitting CTRL-C, thereby potentially resulting in corrupted
	# data being written. Permit these operations to finish before exiting.
	signal.signal(signal.SIGINT, sigterm_handler)

	run_thread = threading.Thread(target=run, args=(safe_to_exit,))
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

if __name__ == "__main__":
	main()
