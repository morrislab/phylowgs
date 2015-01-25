This Python/C++ code is the accompanying software for the paper:
Amit G. Deshwar, Shankar Vembu, Christina K. Yung, Gun Ho Jang, Lincoln Stein, Quaid Morris,
Reconstructing subclonal composition and evolution from whole genome sequencing of tumors.


######################################################################

DEPENDENCIES:
SciPy - http://www.scipy.org/
ETE2 - http://ete.cgenomics.org/
GSL- http://www.gnu.org/software/gsl/

#######################################################################

PREPARING THE INPUT FILE:
Input file format:
The input to evolve.py should be two tab-delimited text files, one for SSM data and one for CNV data. The required column headers are:
SSM DATA
-- id: identifier for each SSM (each row should have a unique identifier)
-- gene: name for each somatic variant
-- a: number of reference allele read counts on the variant locus 
-- d: total number of reads at the locus 
-- mu_r: fraction of expected reference allele sampling from reference population (e.g. if it is an A->T somatic mutation at the locus, the genotype of the reference population should be AA, so the mu_r should be 1-sequencing error rate)
-- mu_v: fraction of expected reference allele sampling from variant population (e.g. if it is an A->T somatic mutation at the locus, copy number is 2 and the expected genotype is AT for the variant population, then the expected fraction of expected reference should be 0.5)

CNV DATA
-- cnv: identifier for each CNV (each row should have a unique identifier)
-- a: number of reference allele read counts on the variant locus 
-- d: total number of reads at the locus 
-- ssms: ssms that overlap with this cnv, each entry is a triplet consisting of ssm id, maternal and paternal copy number, separated by semicolon.



#######################################################################
COMPILING THE C++ FILES:
To compile the C++ file, execute the following command:
g++ -o mh.o  mh.cpp  util.cpp `gsl-config --cflags --libs`

RUNNING THE SOFTWARE:
1. To run the code on the sample data set "ssm_data.txt" and "cnv_data.txt" execute the following command:

python evolve.py 'ssm_data.txt' 'cnv_data.txt' 'trees' 'top_k_trees' 'clonal_frequencies' 'llh_trace' 100 5000 1

# 'ssm_data.txt': file name of the input ssm data
# 'cnv_data.txt': file name of the input cnv data
# 'trees': output folder name where the MCMC trees/samples are saved 
# 'top_k_trees': output file name to save top-k trees in text format. 
# 'clonalFrequencies': output file to save clonal frequencies
# 'llh_trace': output file name to save log likelihood trace
The last three arguments are the number of MCMC samples, number of Metropolis-Hastings iterations, and a random seed number for initializing the MCMC sampler.

2. To generate the posterior trees in PDF/latex format, execute the following commands. The latex files are saved in folder 'latex'

 python posterior_trees.py 'trees' 'ssm_data.txt' 'cnv_data.txt'
 rm *.aux
 rm *.log
 mv *.pdf ./latex/
  
# 'trees': the folder name where the MCMC trees/samples are located (saved from the above command)
# 'ssm_data.txt': file name of the input ssm data
# 'cnv_data.txt': file name of the input cnv data
#######################################################################


#######################################################################
LICENSE:

Copyright (C) 2015 Quaid Morris

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
