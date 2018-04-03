import argparse
import os
import subprocess
import random

def parse_args():
    parser = argparse.ArgumentParser(
		description='Run multiple chains of PhyloWGS and choose run with best tree',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
    );
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
    parser.add_argument('-r', '--random-seeds', dest='random_seeds', default=[], type=list,
		help='Random seeds for initializing MCMC sampler')
    parser.add_argument('-n', '--num-chains', dest='num_chains', default=10, type=int,
		help='Number of chains to run concurrently')
    #parser.add_argument('-t', '--tmp-dir', dest='tmp_dir',
	#	help='Path to directory for temporary files')
    #parser.add_argument('-p', '--params', dest='params_file',
	#	help='JSON file listing run parameters, generated by the parser')
    parser.add_argument('ssm_file',
		help='File listing SSMs (simple somatic mutations, i.e., single nucleotide variants. For proper format, see README.md.')
    parser.add_argument('cnv_file',
		help='File listing CNVs (copy number variations). For proper format, see README.md.')
    args = parser.parse_args()
    if not args.random_seeds:
        random.seed(0);
        args.random_seeds = [random.randint(1,2**32) for i in range(args.num_chains)];
    
    return args

def check_args(args):
    if len(args.random_seeds) != args.num_chains:
        raise ValueError("Number of chains is not equal to the number of input seeds.");

def create_directory(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname);   

def run(args,chain_index,working_dir,output_dir):
    bashCommand = working_dir + "/evolve.py" + \
        " -b " + str(args.write_backups_every) + \
        " -S " + str(args.write_state_every) + \
        " -k " + str(args.top_k_trees) + \
        " -f " + str(args.clonal_freqs) + \
        " -B " + str(args.burnin_samples) + \
        " -s " + str(args.mcmc_samples) + \
        " -i " + str(args.mh_iterations) + \
        " -r " + str(args.random_seeds[chain_index]) + \
        " " + working_dir + "/" + args.ssm_file + \
        " " + working_dir + "/" + args.cnv_file;
    print "Starting chain " + str(chain_index);
    #Note, piping the output may be useful later if we wish to append a "run#:" to the output
    #process = subprocess.Popen(bashCommand,stdout=subprocess.PIPE, shell=True, cwd=output_dir)
    process = subprocess.Popen(bashCommand, shell=True, cwd=output_dir)
    return process

def run_chains(args):
    working_dir = os.path.dirname(os.path.realpath(__file__));
    create_directory(working_dir+"/multevolve_runs");
    run_subprocesses = [];
    for chain_index in range(args.num_chains):
        output_dir = working_dir+"/multevolve_runs/run_"+str(chain_index);
        create_directory(output_dir);
        process = run(args,chain_index,working_dir,output_dir);
        run_subprocesses.append(process);
    
    while True:
        all_dead = all([run_subprocesses[i].poll()!=None for i in range(args.num_chains)]);
        if all_dead:
            break;
            
def main():
    #To do:
    # - Once all runs are done, examine data, keep run with highest likelihood, delete the others. Do this later, will want to look at all runs while testing
    # - Comment all functions.
    # - improve outputs. Should be able to catch outputs from each subprocess and append a "run#" or "chain#".
    # - had to remove some input parameters (params file and temp_dir) as there is no default value in the parser, but rather is handled separately in evolve.
    #   - Can add an if statement in here to catch if these values are none, and omit them from the bash command.
    args = parse_args();
    check_args(args);
    run_chains(args);

if __name__ == "__main__":
    main()