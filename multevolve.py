import argparse
import os
import subprocess
import random
import signal
import zipfile
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(
      description='Concurrently run multiple chains of PhyloWGS.',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter
      )
    parser.add_argument('-n', '--num-chains', dest='num_chains', default=10, type=int,
		  help='Number of chains to run concurrently')
    parser.add_argument('-B', '--burnin-samples', dest='burnin_samples', default=1000, type=int,
		  help='Number of burnin samples')
    parser.add_argument('-s', '--mcmc-samples', dest='mcmc_samples', default=2500, type=int,
		  help='Number of MCMC samples')
    parser.add_argument('-r', '--random-seeds', dest='random_seeds', default=[], type=list,
		  help='Random seeds for initializing MCMC')
    parser.add_argument('-od', '--output-directory', dest='output_directory', default='', type=str,
		  help='Directory where results from each chain will be saved. If directory does not exist, will attempt to create it here. (Default = "current_directory/multevolve_chains")')
    parser.add_argument('-sf','--ssm-file',dest='ssm_file',
		help='File listing SSMs (simple somatic mutations, i.e., single nucleotide variants.')
    parser.add_argument('-cf','--cnv-file',dest='cnv_file',
		help='File listing CNVs (copy number variations).')
    args, evolve_args = parser.parse_known_args()
    args = dict(args._get_kwargs());

    return args, evolve_args

def check_args(args):
    #Set default values for arguments that require function calls to calculate
    if not args['random_seeds']:
        random.seed(0);
        args['random_seeds'] = [random.randint(1,2**32) for i in range(args['num_chains'])];
    if not args['output_directory']:
        args['output_directory'] = os.path.join(os.getcwd(),"multevolve_chains");
        create_directory(args['output_directory']);

    #Make sure the arguments make sense. Right now just have to check that the 
    #list of random seeds, if this was provided by the user, has length = num_chains.
    if len(args['random_seeds']) != args['num_chains']:
        raise ValueError("Number of chains is not equal to the number of input seeds.");
    
    return args;

def create_directory(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname);   

def run_chains(args,evolve_args):
    '''
    Determine location of evolve.py (same directory as this script), location of the ssm
    and cnv files, and create the output directories for each chain. Create a subprocess
    for each chain so that they may all run at the same time.
    '''
    current_dir = os.getcwd(); 
    working_dir = os.path.dirname(os.path.realpath(__file__));
    processes = [];
    out_dirs = [];
    for chain_index in range(args['num_chains']):
        output_dir = os.path.join(args['output_directory'],"chain_"+str(chain_index));
        out_dirs.append(output_dir);
        create_directory(output_dir);
        process = run(args,evolve_args,chain_index,working_dir,current_dir,output_dir);
        processes.append(process);
    watch_chains(args,processes);
    return out_dirs

def run(args,evolve_args,chain_index,working_dir,current_dir,output_dir):
    '''
    Start a new subprocess for every call to evolve. Return the subprocess 
    so that we can capture its outputs and see if it is complete.
    '''
    bashCommand = ["python2", \
        os.path.join(working_dir, "evolve.py"), \
        "-B", str(args['burnin_samples']), \
        "-s", str(args['mcmc_samples']), \
        "-r", str(args['random_seeds'][chain_index]),\
        os.path.join(current_dir,args['ssm_file']), \
        os.path.join(current_dir,args['cnv_file'])];
    bashCommand = bashCommand + list(evolve_args);
    print "Starting chain " + str(chain_index);
    process = subprocess.Popen(bashCommand,stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=output_dir, universal_newlines=True,bufsize=1)
    return process

def watch_chains(args,processes):
    '''
    This watches the subprocesses as they run. Will check to see if they
    are complete and if so, move on to the next step; and if not, then will
    capture their outputs and inform the user of each processes progression.
    '''
    def stdout_read_alarm_handler(signum, frame):
        raise Exception("subprocess.stdout.readline() took too long to complete.")
    signal.signal(signal.SIGALRM,stdout_read_alarm_handler);

    num_chains = args['num_chains']
    progression_text = ['\n']*len(processes);
    print "".join(progression_text)
    while True:
        #Check to see if all processes are done running and if so, exit.
        all_dead = all([processes[i].poll()!=None for i in range(num_chains)]);
        if all_dead:
            break;
        
        # Capture the output from each process, modify it, and output it. 
        # Note: process.stdout returns a file handle from which you can read outputs, however if you
        # call readline and there is no new information to read, it will wait until there is until
        # returning, which is really annoying. So I set up a timer here to throw an exception if
        # more than 0.1 seconds passes so that we can move on to read info from another process that may already 
        # have output to read.
        other_text = [];
        for chain_index in range(num_chains):
            #If the process is finished, then will return an empty string when calling readline(), skip this.
            if processes[chain_index].poll()!=None:
                continue
            try:
              signal.setitimer(signal.ITIMER_REAL,0.1)
              new_line = processes[chain_index].stdout.readline()
              signal.setitimer(signal.ITIMER_REAL,0)
            except Exception, e:
              signal.setitimer(signal.ITIMER_REAL,0)
              continue
            if (len(new_line.split()) > 2) and new_line.split()[2].lstrip('-').isdigit():
              trees_done = int(new_line.split()[2])+args['burnin_samples']
              total_trees = args['burnin_samples']+args['mcmc_samples']
              percent_complete = 100*trees_done/total_trees
              progression_text[chain_index] = "chain{}: {}/{} - {}% complete\n".format(chain_index, trees_done, total_trees, percent_complete)
            else: 
              other_text.append("chain{}: {}".format(chain_index,new_line))
        print "\033[F"*(num_chains+1), # Move cursor up to line that starts telling us about chain progression. Want to overwrite those lines.
        print "".join(other_text),
        print " "*70
        print "".join(progression_text),

def determine_each_chains_highest_likelihood(chain_dirs):
    '''
    Examines all of the trees output by each chain and reports
    the highest likelihood of all of the trees.
    '''
    best_likelihoods = [];
    for chain_dir in chain_dirs:
        this_best_likelihood = -(2**32);
        tree_zip_file = zipfile.ZipFile(os.path.join(chain_dir,'trees.zip'), mode = 'r')
        for tree_name in tree_zip_file.namelist():
            if "tree" in tree_name:
                likelihood = float(tree_name.split('_')[-1]) #likelihood is in the names of the trees, just use that.
                if likelihood > this_best_likelihood:
                    this_best_likelihood = likelihood;
        best_likelihoods.append(this_best_likelihood);
    return best_likelihoods

def merge_best_chains(args,chain_dirs,best_likelihoods):
    '''
    Determines which chains are the best and merges them together into one trees.zip
    file that can be input into write_results.
    A chain counts, for now, as being one of the best if the highest likelihood of all 
    of it's trees is within 10% of the highest likelihood of all of the trees calculated
    across all chains.
    '''
    overall_best_likelihood = max(best_likelihoods);
    #Let's say, arbitrarily for now, that any best likelihoods that are within 10% of the best best, have their chains merged into one.
    chains_to_include = [ind for ind,bl in enumerate(best_likelihoods) if abs((bl-overall_best_likelihood)/overall_best_likelihood) < 0.1];
    out_dir = os.path.join(args['output_directory'],'merged_best_chains')
    create_directory(out_dir)
    combined_tree_zipfile = zipfile.ZipFile(os.path.join(out_dir,"trees.zip"), mode='w', compression=zipfile.ZIP_DEFLATED, allowZip64=True);
    print "Merging best chains:"
    tree_index = 0;
    for chain_idx in chains_to_include:
        print "  merging chain " + str(chain_idx) + "...";
        chain_dir = chain_dirs[chain_idx];
        this_zip = zipfile.ZipFile(os.path.join(chain_dir,"trees.zip"), mode='r');
        this_zips_files = this_zip.namelist();
        files_to_include = [(filename, this_zip.read(filename)) for filename in this_zips_files if "tree" in filename]
        for file in files_to_include:
            filename_components = file[0].split("_")
            filename_components[1] = str(tree_index)
            filename = "_".join(filename_components)
            combined_tree_zipfile.writestr(filename, file[1])
            tree_index += 1

def main():
    args,evolve_args = parse_args();
    check_args(args);
    chain_dirs = run_chains(args,evolve_args);
    each_chains_best_likelihoods = determine_each_chains_highest_likelihood(chain_dirs);
    merge_best_chains(args, chain_dirs, each_chains_best_likelihoods);

if __name__ == "__main__":
    main()