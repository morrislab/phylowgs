import argparse
import os
import subprocess
import random
import sys
import signal

def parse_args():
    parser = argparse.ArgumentParser(
      description='Run multiple chains of PhyloWGS.',
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
    args, evolve_args = parser.parse_known_args()
    args = dict(args._get_kwargs())
    if not args['random_seeds']:
        random.seed(0);
        args['random_seeds'] = [random.randint(1,2**32) for i in range(args['num_chains'])];
    return args, evolve_args

def check_args(args):
    #Just making sure the arguments input make sense. Right now just have to check that 
    #the list of random seeds input, if this was provided, has length = num_chains.
    if len(args['random_seeds']) != args['num_chains']:
        raise ValueError("Number of chains is not equal to the number of input seeds.");

def create_directory(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname);   

def run(args,evolve_args,chain_index,working_dir,current_dir,output_dir):
    '''
    Start a new subprocess for every run of evolve. Return the subprocess 
    so that we can capture it's outputs and see if it's complete.
    '''
    bashCommand = [working_dir + "/evolve.py", \
        " -B ", str(args['burnin_samples']), \
        " -s ", str(args['mcmc_samples']), \
        " -r ", str(args['random_seeds'][chain_index])];
    bashCommand = bashCommand + list(evolve_args);
    print "Starting chain " + str(chain_index);
    process = subprocess.Popen(bashCommand,stdout=subprocess.PIPE, shell=True,stderr=subprocess.STDOUT, cwd=output_dir, universal_newlines=True,bufsize=1)
    return process

def watch_runs(args,processes):
    '''
    This watches the subprocesses as they run. Will check to see if they
    are complete and if so, end this script; and if not, then will capture
    their outputs and inform the user of each processes progression.
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
        # Note: process.stdout returns a file handle from which you can read outputs, however if
        # call readline and there is no new information to read, it will wait until there is until
        # moving on, which is really annoying. So I set up an alarm here to throw and exception if
        # more than 1 second passes so that we can read info from another process that may already 
        # have more output.
        other_text = [];
        for chain_index in range(num_chains):
            try:
              signal.alarm(1)
              new_line = processes[chain_index].stdout.readline()
              signal.alarm(0)
            except Exception, e:
              signal.alarm(0)
              continue
            if new_line.split()[2].lstrip('-').replace('.','',1).isdigit():
              trees_done = int(new_line.split()[2])+args['burnin_samples']
              total_trees = args['burnin_samples']+args['mcmc_samples']
              percent_complete = (100*trees_done)/total_trees
              progression_text[chain_index] = "chain{}: {}/{} - {}% complete\n".format(chain_index, trees_done, total_trees, percent_complete)
            else: 
              other_text.append("chain{}: {}".format(chain_index,new_line))
        print "\033[F"*(num_chains+1), # Cursor up to line that starts telling us about run progression. Want to overwrite those lines.
        print "".join(other_text),
        print " "*50
        print "".join(progression_text),
        

def run_chains(args,evolve_args):
    current_dir = os.getcwd();
    working_dir = os.path.dirname(os.path.realpath(__file__));
    create_directory(current_dir+"/multevolve_runs");
    processes = [];
    for chain_index in range(args['num_chains']):
        output_dir = current_dir+"/multevolve_runs/run_"+str(chain_index);  
        create_directory(output_dir);
        process = run(args,evolve_args,chain_index,working_dir,current_dir,output_dir);
        processes.append(process);
    watch_runs(args,processes);
    
            
def main():
    #To do:
    # - Once all runs are done, examine data, keep run with highest likelihood, delete the others. Do this later, will want to look at all runs while testing
    #    - methinks this will be a separate script.
    
    args,evolve_args = parse_args();
    check_args(args);
    run_chains(args,evolve_args);

if __name__ == "__main__":
    main()