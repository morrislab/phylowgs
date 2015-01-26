This Python/C++ code is the accompanying software for the paper:
Amit G. Deshwar, Shankar Vembu, Christina K. Yung, Gun Ho Jang, Lincoln Stein, Quaid Morris,
Reconstructing subclonal composition and evolution from whole genome sequencing of tumors.


######################################################################


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
USAGE:

  1. Install dependencies.

    # Install Python 2 versions of NumPy (www.numpy.org) and SciPy (www.scipy.org).
    # Install Python 2 version of ETE2 (e.g.: pip2 install --user ete2).
    # Install GSL (http://www.gnu.org/software/gsl/).

  2. Compile the C++ file.

    g++ -o mh.o  mh.cpp  util.cpp `gsl-config --cflags --libs`

  3. Run PhyloWGS.

    # Minimum invocation on sample data set: python2 evolve.py ssm_data.txt cnv_data.txt

    # All options:

        usage: evolve.py [-h] [-t TREES] [-k TOP_K_TREES] [-f CLONAL_FREQS]
                         [-l LLH_TRACE] [-s MCMC_SAMPLES] [-i MH_ITERATIONS]
                         [-r RANDOM_SEED]
                         ssm_file cnv_file

        positional arguments:
          ssm_file              File listing SSMs (simple somatic mutations, i.e.,
                                single nucleotide variants. For proper format, see
                                README.txt.
          cnv_file              File listing CNVs (copy number variations). For proper
                                format, see README.txt.

        optional arguments:
          -h, --help            show this help message and exit
          -t TREES, --trees TREES
                                Output folder name where the MCMC trees/samples are
                                saved (default: trees)
          -k TOP_K_TREES, --top-k-trees TOP_K_TREES
                                Output file name to save top-k trees in text format
                                (default: top_k_trees)
          -f CLONAL_FREQS, --clonal-freqs CLONAL_FREQS
                                Output file to save clonal frequencies (default:
                                clonalFrequencies)
          -l LLH_TRACE, --llh-trace LLH_TRACE
                                Output file name to save log likelihood trace
                                (default: llh_trace)
          -s MCMC_SAMPLES, --mcmc-samples MCMC_SAMPLES
                                Number of MCMC samples (default: 2500)
          -i MH_ITERATIONS, --mh-iterations MH_ITERATIONS
                                Number of Metropolis-Hastings iterations (default:
                                5000)
          -r RANDOM_SEED, --random-seed RANDOM_SEED
                                Random seed for initializing MCMC sampler. If
                                unspecified, choose random seed automatically.
                                (default: None)


  4. Generate the posterior trees in PDF/latex format. The LaTeX files are
     saved in folder 'latex'.

       python posterior_trees.py 'trees' 'ssm_data.txt' 'cnv_data.txt'
       rm *.aux
       rm *.log
       mv *.pdf ./latex/
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
