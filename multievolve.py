import argparse
import os
import subprocess
import random
import zipfile
import numpy as np
import sys
import re
from util2 import logmsg
import Queue
import threading
from collections import defaultdict

def create_directory(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def logsumexp(a):
    #numpy has a logaddexp() but it only takes two values and that's stupid so
    #I'm just going to create my own logsumexp here.
    max_a = np.max(a)
    result = max_a + np.log(np.sum([np.exp(i-max_a) for i in a]))
    return result

def parse_args():
    parser = argparse.ArgumentParser(
      description='Concurrently run multiple chains of PhyloWGS.',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n', '--num-chains', dest='num_chains', default=10, type=int,
          help='Number of chains to run concurrently')
    parser.add_argument('-r', '--random-seeds', dest='random_seeds', default=[], type=list,
          help='Random seeds for initializing MCMC')
    parser.add_argument('-if', '--chain-inclusion-factor', dest='chain_inclusion_factor', default=1.5, type=float,
          help='Factor for determining which chains will be included in the output "merged" folder. ' \
               'Default is 1.5, meaning that the sum of the likelihoods of the trees found in each chain must ' \
               'be greater than 1.5x the maximum of that value across chains. Setting this value = inf ' \
               'includes all chains and setting it = 1 will include only the best chain.')
    parser.add_argument('-od', '--output-directory', dest='output_directory', default='', type=str,
          help='Directory where results from each chain will be saved. If directory does not exist, ' \
               'will attempt to create it here. (Default = "working_directory/multievolve_chains")')
    # Send unrecognized arguments to evolve.py.
    known_args, other_args = parser.parse_known_args()
    known_args = dict(known_args._get_kwargs())
    return known_args, other_args

def check_args(args):
    #Set default values for arguments that require function calls to calculate
    if not args['random_seeds']:
        random.seed(0)
        args['random_seeds'] = [random.randint(1,2**32) for i in range(args['num_chains'])]
    if not args['output_directory']:
        args['output_directory'] = os.path.join(os.getcwd(),"multievolve_chains")
        create_directory(args['output_directory'])

    #Make sure the arguments make sense. Right now just have to check that the
    #list of random seeds, if this was provided by the user, has length = num_chains.
    if len(args['random_seeds']) != args['num_chains']:
        raise ValueError("Must specify random seeds for every chain")
    return args

def run_chains(args,evolve_args):
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
        output_dir = os.path.join(args['output_directory'],"chain_"+str(chain_index))
        out_dirs.append(output_dir)
        create_directory(output_dir)
        process = run(args,evolve_args,chain_index,app_dir,working_dir,output_dir)
        processes.append(process)
    watch_chains(processes)
    return out_dirs

def run(args,evolve_args,chain_index,app_dir,working_dir,output_dir):
    '''
    Start a new subprocess for every call to evolve. Return the subprocess
    so that we can capture its outputs and see if it is complete.
    '''
    cmd = [
        sys.executable,
        os.path.join(app_dir, "evolve.py"),
        '--output-dir', output_dir,
    ]
    cmd = cmd + list(evolve_args)
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
    delay = 0.5

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

            try:
                line = Q.get(timeout=delay)
            except Queue.Empty:
                continue
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
        logSumLHs.append(logsumexp(logLHs))

    # Check below assumes that LLH < 0, which it should always be.
    logSumLHs = np.array(logSumLHs)
    assert np.all(logSumLHs < 0)
    bestLogSumLH = np.max(logSumLHs)
    chains_to_merge = [i for i,logSumLH in enumerate(logSumLHs) if logSumLH > (chain_inclusion_factor*bestLogSumLH)]
    return chains_to_merge

def merge_best_chains(args,chain_dirs,chains_to_merge):
    '''
    Determines which chains are the best and merges them together into one trees.zip
    file that can be input into write_results.
    A chain counts, for now, as being one of the best if the highest likelihood of all
    of it's trees is within 10% of the highest likelihood of all of the trees calculated
    across all chains.
    '''
    out_dir = os.path.join(args['output_directory'],'merged_best_chains')
    create_directory(out_dir)
    if os.path.isfile(os.path.join(out_dir,"trees.zip")):
        logmsg("Merged trees.zip file already exists. To create a new merged trees.zip, remove the existing one first.")
        return
    combined_tree_zipfile = zipfile.ZipFile(os.path.join(out_dir,"trees.zip"), mode='w', compression=zipfile.ZIP_DEFLATED, allowZip64=True)
    logmsg("Merging best chains:")
    tree_index = 0
    for chain_idx in chains_to_merge:
        logmsg("  merging chain {} ...".format(chain_idx))
        chain_dir = chain_dirs[chain_idx]
        this_zip = zipfile.ZipFile(os.path.join(chain_dir,"trees.zip"), mode='r')
        this_zips_files = this_zip.namelist()
        files_to_include = [(filename, this_zip.read(filename)) for filename in this_zips_files if filename.startswith("tree")]
        for file in files_to_include:
            #First we need to reindex the tree
            filename_components = file[0].split("_")
            filename_components[1] = str(tree_index)
            filename = "_".join(filename_components)
            combined_tree_zipfile.writestr(filename, file[1])
            tree_index += 1
    #Don't forget the "params.json" and "cnv_logical_physical_mapping.json" files. They should all be the same in each
    #chains zip file. So just take the last one used and insert it.
    combined_tree_zipfile.writestr("cnv_logical_physical_mapping.json", this_zip.read("cnv_logical_physical_mapping.json"))
    combined_tree_zipfile.writestr("params.json", this_zip.read("params.json"))
    logmsg("Chain merging complete.")
    logmsg("You can remove unneeded intermediate files via the shell command `rm /path/to/output/dir/multievolve_chains/chain_*/trees.zip`")

def main():
    args,evolve_args = parse_args()
    check_args(args)
    chain_dirs = run_chains(args,evolve_args)
    chains_to_merge = determine_chains_to_merge(chain_dirs, args['chain_inclusion_factor'])
    merge_best_chains(args, chain_dirs, chains_to_merge)

if __name__ == "__main__":
    main()
