import argparse
import os
import subprocess
import zipfile
import numpy as np
import sys
import re
from util2 import logmsg
import Queue
import threading
import time
import scipy.misc
import hashlib
from collections import defaultdict

def create_directory(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def parse_args():
    parser = argparse.ArgumentParser(
      description='Concurrently run multiple MCMC chains of PhyloWGS. ' +
      'All options that evolve.py accepts may also be specified here. To list those arguments, run `%s evolve.py --help`.' % sys.executable,
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n', '--num-chains', dest='num_chains', default=4, type=int,
          help='Number of chains to run concurrently')
    parser.add_argument('-r', '--random-seeds', dest='random_seeds', type=int, nargs='+',
          help='Space-separated random seeds with which to initialize each chain. Specify one for each chain.')
    parser.add_argument('-I', '--chain-inclusion-factor', dest='chain_inclusion_factor', default=1.5, type=float,
          help='Factor for determining which chains will be included in the output "merged" folder. ' \
               'Default is 1.5, meaning that the sum of the likelihoods of the trees found in each chain must ' \
               'be greater than 1.5x the maximum of that value across chains. Setting this value = inf ' \
               'includes all chains and setting it = 1 will include only the best chain.')
    parser.add_argument('-O', '--output-dir', dest='output_dir', default='chains',
          help='Directory where results from each chain will be saved. We will create it if it does not exist.')
    # Ideally, I wouldn't specify `ssm_file` or `cnv_file` as multievolve.py
    # arguments, since I don't need them here -- I just want to pass them
    # through to evolve.py. But then printing the help is confusing, as you
    # don't realize that you should pass them as arguments to multievolve. So,
    # I should specify them here -- otherwise, you can invoke multievolve.py
    # with no ssm_data.txt or cnv_data.txt, and it will dutifully invoke
    # multiple copies of evolve.py with no input files (with the evolve.py runs
    # immediately failing.
    #
    # Since --ssms and --cnvs are required, it would be better to make them
    # positional arguments. But this screws up parsing with parse_known_args()
    # -- it mixes up unknown and known arguments if I make ssm_file and cnv_file positional, then call
    # `python2 ../multievolve.py -n2 -s 5 -B 3 ../ssm_data.txt ../cnv_data.txt`.
    # To fix this, just make all known arguments for multievolve (required) optional arguments.
    #
    # Much of this mess arises from having a wrapper script to start multiple
    # chains. In a better world, evolve.py would natively support multiple
    # chains, with the user able to choose to have only a single chain if she
    # desires.
    parser.add_argument('--ssms', dest='ssm_file', required=True,
            help='File listing SSMs (simple somatic mutations, i.e., single nucleotide variants. For proper format, see README.md.')
    parser.add_argument('--cnvs',dest='cnv_file', required=True,
            help='File listing CNVs (copy number variations). For proper format, see README.md.')

    # Send unrecognized arguments to evolve.py.
    known_args, other_args = parser.parse_known_args()
    return dict(known_args._get_kwargs()), other_args

def check_args(args):
    args['output_dir'] = os.path.abspath(args['output_dir'])
    create_directory(args['output_dir'])

    #Make sure the arguments make sense. Right now just have to check that the
    #list of random seeds, if this was provided by the user, has length = num_chains.
    if args['random_seeds'] is not None and len(args['random_seeds']) != args['num_chains']:
        raise ValueError("Must specify exactly one seed for each of %s chain(s). You specified %s seed(s)." % (
            args['num_chains'], len(args['random_seeds'])))
    return args

def run_chains(args, other_args):
    '''
    Determine location of evolve.py (same directory as this script), location of the ssm
    and cnv files, and create the output directories for each chain. Create a subprocess
    for each chain so that they may all run at the same time.
    '''
    working_dir = os.getcwd()
    app_dir = os.path.dirname(os.path.realpath(__file__))
    processes = []
    out_dirs = []
    for chain_index in range(args['num_chains']):
        output_dir = os.path.join(args['output_dir'],"chain_"+str(chain_index))
        out_dirs.append(output_dir)
        create_directory(output_dir)
        process = run_chain(chain_index, args['ssm_file'], args['cnv_file'], app_dir, working_dir, output_dir, args['random_seeds'], other_args)
        processes.append(process)
    watch_chains(processes)
    return out_dirs

def run_chain(chain_index, ssm_fn, cna_fn, app_dir, working_dir, output_dir, seeds, other_args):
    '''
    Start a new subprocess for every call to evolve. Return the subprocess
    so that we can capture its outputs and see if it is complete.
    '''
    cmd = [
        sys.executable,
        os.path.join(app_dir, "evolve.py"),
        '--output-dir', output_dir,
    ]
    if seeds is not None:
        cmd += ['--random-seed', str(seeds[chain_index])]
    cmd += [
        ssm_fn,
        cna_fn,
    ]
    cmd = cmd + list(other_args)

    logmsg("Starting chain %s" % chain_index)
    # bufsize=1 and universal_newlines=True open stdout in line-buffered text
    # mode, rather than binary stream.
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
        universal_newlines=True,
        close_fds=True,
        cwd=working_dir,
    )
    return process

def parse_status(line):
    status = {}
    fields = line.split(' ')
    for F in fields:
        K, V = F.split('=', 1)
        status[K] = V
    return status

def enqueue_output(out, queue):
    for line in iter(out.readline, b''):
        queue.put(line)
    out.close()

def make_queues(processes):
    queues = []

    for P in processes:
        Q = Queue.Queue()
        T = threading.Thread(target=enqueue_output, args=(P.stdout, Q))
        T.daemon = True # Thread dies with program
        T.start()
        queues.append(Q)

    return queues

def watch_chains(processes):
    # Based on https://stackoverflow.com/a/4896288. TL;DR: non-blocking reads
    # from stdout on subprocesses are really painful. This is the cleanest and
    # most reliable mechanism I've come across for resolving them.
    num_chains = len(processes)
    last_lines_were_status = False
    status = {idx: {'status': 'initializing'} for idx in range(num_chains)}
    chain_stdout = defaultdict(list)
    delay = 0.05

    queues = make_queues(processes)
    while True:
        # All are done.
        if set([S['status'] for S in status.values()]) == set(['done']):
            break
        for chain_index, P, Q in zip(range(num_chains), processes, queues):
            if status[chain_index]['status'] == 'done':
                continue
            exit_code = processes[chain_index].poll()
            if exit_code is not None:
                # Note we still finish the rest of this loop iteration, which
                # lets us print the process' final output.
                status[chain_index] = {'status': 'done', 'exit_code': exit_code}

            while True:
                # Use loop so that we retrieve as many lines as are available
                # from the process.
                try:
                    line = Q.get(timeout=delay)
                except Queue.Empty:
                    break
                else:
                    chain_stdout[chain_index].append(line.strip())

        for chain_index in sorted(chain_stdout.keys()):
            for line in chain_stdout[chain_index]:
                # Strip existing timestamp.
                if re.match(r'^\[\d{4}-', line):
                    line = line[line.index(']')+1:].strip()
                if line.startswith('iteration='):
                    status_line = parse_status(line)
                    status[chain_index] = parse_status(line)
                    status[chain_index]['status'] = 'running'
                    status[chain_index]['percent_complete'] = '{:.2f}%'.format(100 * float(status[chain_index]['trees_sampled']) / float(status[chain_index]['total_trees']))
                else:
                    if len(line) > 0:
                        logmsg("chain={} {}".format(chain_index, line))
                        last_lines_were_status = False
        chain_stdout = defaultdict(list)

        if last_lines_were_status and sys.stdout.isatty():
            print("\033[2K\033[1A" * (num_chains + 1)) # Move cursor up to line that starts telling us about chain progression. Want to overwrite those lines.
        for cidx in sorted(status.keys()):
            if status[cidx]['status'] == 'running':
                keys = ('trees_sampled', 'total_trees', 'percent_complete')
            elif status[cidx]['status'] == 'done':
                keys = ('exit_code',)
            else:
                keys = tuple()
            status_msg = ' '.join(['{}={}'.format(K, status[cidx][K]) for K in ('status',) + keys])
            logmsg('chain={} {}'.format(cidx, status_msg))
            last_lines_were_status = True
        time.sleep(1)

def determine_chains_to_merge(chain_dirs,chain_inclusion_factor):
    '''
    Examines all of the trees output by each chain and reports which chains should
    be merged. Chains will meet the criteria if the log(sum(all_tree_likelihoods))
    is within some factor of the maximum of that value across chains.
    '''
    logSumLHs = []
    for chain_dir in chain_dirs:
        logLHs = []
        tree_zip_file = zipfile.ZipFile(os.path.join(chain_dir,'trees.zip'), mode = 'r')
        for tree_name in tree_zip_file.namelist():
            if tree_name.startswith("tree"):
                #the logged likelihood is in the names of the trees, just use that.
                logLHs.append(float(tree_name.split('_')[-1]))
        logSumLHs.append(scipy.misc.logsumexp(logLHs))

    # Check below assumes that LLH < 0, which it should always be. We need this
    # assumption for the idea that a "slightly worse" chain has a "slightly
    # more negative" LH to work.
    logSumLHs = np.array(logSumLHs)
    assert np.all(logSumLHs < 0)
    bestLogSumLH = np.max(logSumLHs)

    included_chains = []
    excluded_chains = []
    for cidx, logsumlh in enumerate(logSumLHs):
        if logsumlh >= chain_inclusion_factor * bestLogSumLH:
            included_chains.append((cidx, logsumlh))
        else:
            excluded_chains.append((cidx, logsumlh))

    assert len(included_chains) >= 1
    return (included_chains, excluded_chains)

def merge_best_chains(out_dir, chain_dirs, included_chains, excluded_chains):
    '''
    Determines which chains are the best and merges them together into one trees.zip
    file that can be input into write_results.
    A chain counts, for now, as being one of the best if the highest likelihood of all
    of it's trees is within 10% of the highest likelihood of all of the trees calculated
    across all chains.
    '''
    combined_fn = os.path.join(out_dir,"trees.zip")
    if os.path.isfile(combined_fn):
        logmsg("Merged trees.zip file already exists. To create a new merged trees.zip, remove the existing one first.")
        return

    combined_tree_zipfile = zipfile.ZipFile(combined_fn, mode='w', compression=zipfile.ZIP_DEFLATED, allowZip64=True)
    logmsg("Including chains {}".format(' '.join(['{}={}'.format(cidx, logsumlh) for cidx, logsumlh in included_chains])))
    if len(excluded_chains) > 0:
        logmsg("Excluding chains {}".format(' '.join(['{}={}'.format(cidx, logsumlh) for cidx, logsumlh in excluded_chains])))
    else:
        logmsg('Not excluding any chains')
    tree_index = 0
    zip_paths = []
    others = defaultdict(dict)

    for chain_idx, _ in included_chains:
        chain_dir = chain_dirs[chain_idx]
        zip_path = os.path.abspath(os.path.join(chain_dir, 'trees.zip'))
        zip_paths.append(zip_path)
        this_zip = zipfile.ZipFile(zip_path, mode='r')
        this_zips_files = this_zip.namelist()

        is_tree_file = lambda fn: fn.startswith('tree')
        should_include_other = lambda fn: not fn.startswith('burnin')
        tree_fns  = [fn for fn in this_zips_files if is_tree_file(fn)]
        other_fns = [fn for fn in this_zips_files if not is_tree_file(fn) and should_include_other(fn)]

        for fn in tree_fns:
            F = this_zip.read(fn)
            filename_components = fn.split("_")
            # First we need to reindex the tree
            filename_components[1] = str(tree_index)
            fn = "_".join(filename_components)
            combined_tree_zipfile.writestr(fn, F)
            tree_index += 1

        for fn in other_fns:
            others[fn][chain_idx] = hashlib.sha256(this_zip.read(fn)).hexdigest()

    # Assume that any non-tree files are identical across individual trees.zip
    # files. But let's check this assumption.
    for fn in others.keys():
        assert len(set(others[fn].values())) == 1
        combined_tree_zipfile.writestr(fn, this_zip.read(fn))

    write_results_path = os.path.normpath(os.path.join(os.path.dirname(__file__), 'write_results.py'))
    logmsg("Chain merging complete.")
    logmsg("You may remove the following unneeded intermediate files: {}".format(' '.join(zip_paths)))
    logmsg('To write JSON results, please run `{} {} {} {} {} {} {}`'.format(
        sys.executable,
        write_results_path,
        'run_name',
        combined_fn,
        'run_name.summ.json.gz',
        'run_name.muts.json.gz',
        'run_name.mutass.zip'
    ))

def main():
    args,evolve_args = parse_args()
    check_args(args)
    chain_dirs = run_chains(args,evolve_args)
    included_chains, excluded_chains = determine_chains_to_merge(chain_dirs, args['chain_inclusion_factor'])
    merge_best_chains(args['output_dir'], chain_dirs, included_chains, excluded_chains)

if __name__ == "__main__":
    main()
