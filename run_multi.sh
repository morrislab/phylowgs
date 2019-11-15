#!/bin/bash
# run_multi.sh
#SBATCH --partition=exacloud
#SBATCH --output=multievolve-%j.out
#SBATCH --error=multievolve-%j.err
#SBATCH --job-name=run_smchet_multievolve
#SBATCH --gres disk:1024
#SBATCH --mincpus=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=06:00:00

source /home/groups/EllrottLab/activate_conda
ABS_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

sbatch cwltool $ABS_PATH/multievolve.cwl $ABS_PATH/run_multi.json
