#!/bin/bash
# run_parser.sh
#SBATCH --partition=exacloud
#SBATCH --output=parser-%j.out
#SBATCH --error=parser-%j.err
#SBATCH --job-name=run_smchet_parser
#SBATCH --gres disk:1024
#SBATCH --mincpus=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:45:00

source /home/groups/EllrottLab/activate_conda
ABS_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

sbatch cwltool --debug $ABS_PATH/parse_cnvs.cwl $ABS_PATH/parse_cnvs.json
