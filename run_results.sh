#!/bin/bash
# run_results.sh
#SBATCH --partition=exacloud
#SBATCH --output=results-%j.out
#SBATCH --error=results-%j.err
#SBATCH --job-name=run_smchet_results
#SBATCH --gres disk:1024
#SBATCH --mincpus=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:45:00

source /home/groups/EllrottLab/activate_conda
ABS_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

sbatch cwltool --debug $ABS_PATH/write_results.cwl $ABS_PATH/write_results.json
