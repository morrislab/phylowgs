PhyloWGS
========

This Python/C++ code is the accompanying software for the paper [PhyloWGS:
Reconstructing subclonal composition and evolution from whole-genome sequencing
of tumors](http://genomebiology.com/2015/16/1/35), with authors Amit G.
Deshwar, Shankar Vembu, Christina K. Yung, Gun Ho Jang, Lincoln Stein, and Quaid
Morris.


Input files
-----------
The input to evolve.py is two tab-delimited text files -- one for SSM data and
one for CNV data. Please see the files `ssm_data.txt` and `cnv_data.txt`
included with PhyloWGS for examples.

`ssm_data.txt`:

* `id`: identifier for each SSM. Identifiers must start at `s0` and
  increment, so the first data row will have `s0`, the second row `s1`, and so
  forth.
* `gene`: any string identifying the variant -- this need not be a gene name.
  `<chr>_<pos>` (e.g., `2_234577`) works well.
* `a`: number of reference-allele reads at the variant locus.
* `d`: total number of reads at the variant locus.
* `mu_r`: fraction of expected reference allele sampling from the *reference*
  population. E.g., if the tumor has an A->T somatic mutation at the locus,
  the genotype of the reference population should be AA. Thus, `mu_r` should
  be `1 - (sequencing error rate)`. Given the 0.001 error rate in Illumina
  sequencing, setting this column to 0.999 works well.
* `mu_v`: fraction of expected reference allele sampling from *variant*
  population. Suppose an A->T somatic mutation occurred at the locus. `mu_v`
  always uses normal ploidy (i.e., the copy number in non-CNV regions). As
  humans are diploid, copy number will thus always be 2. So, the variant
  population genotype should be AT, meaning we will observe the reference
  allele with frequency `0.5 - (sequencing error rate)`. Given the 0.001
  error rate in Illumina sequencing, setting this column to 0.499 works well.

`cnv_data.txt`:

* `cnv`: identifier for each CNV. Identifiers must start at `c0` and
  increment, so the first data row will have `c0`, the second row `c1`, and so
  forth.
* `a`: number of reference reads covering the CNV.
* `d`: total number of reads covering the CNV. This will be affected by
  factors such as total copy number at the locus, sequencing depth, and the
  size of the chromosomal region spanned by the CNV.
* `ssms`: SSMs that overlap with this CNV. Each entry is a comma-separated
  triplet consisting of SSM ID, maternal copy number, and paternal copy
  number. These triplets are separated by semicolons.


Running PhyloWGS
----------------

1. Install dependencies.

  * Install Python 2 versions of NumPy (www.numpy.org) and SciPy (www.scipy.org).
  * Install Python 2 version of ETE2 (e.g.: `pip2 install --user ete2`).
  * Install GSL (http://www.gnu.org/software/gsl/).

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
                                    Output file where the MCMC trees/samples are saved
                                    (default: trees.zip)
              -k TOP_K_TREES, --top-k-trees TOP_K_TREES
                                    Output file to save top-k trees in text format
                                    (default: top_k_trees)
              -f CLONAL_FREQS, --clonal-freqs CLONAL_FREQS
                                    Output file to save clonal frequencies (default:
                                    clonalFrequencies)
              -l LLH_TRACE, --llh-trace LLH_TRACE
                                    Output file to save log likelihood trace (default:
                                    llh_trace)
              -s MCMC_SAMPLES, --mcmc-samples MCMC_SAMPLES
                                    Number of MCMC samples (default: 2500)
              -i MH_ITERATIONS, --mh-iterations MH_ITERATIONS
                                    Number of Metropolis-Hastings iterations (default:
                                    5000)
              -r RANDOM_SEED, --random-seed RANDOM_SEED
                                    Random seed for initializing MCMC sampler. If
                                    unspecified, choose random seed automatically.
                                    (default: None)

4. Generate the posterior trees in PDF & LaTeX formats. The LaTeX files and
   resulting PDFs are saved in the directory `posterior_trees`.

        python2 posterior_trees.py ssm_data.txt cnv_data.txt


Interpreting output
-------------------

Two primary outputs result from running PhyloWGS:

  * After running `evolve.py`, LaTeX and PDF depictions of the five trees with
    the highest likelihoods are saved in the `top_trees` directory.

  * After running `posterior_trees.py`, LaTeX and PDF depictions of all trees
    are saved in the `posterior_trees` directory, ordered by posterior probability.
    Likelihoods of the trees are ignored -- instead, posterior probabilities
    represent the frequency with which trees bearing the same subclone structure
    and same assignment of mutations to subclones were sampled.

Trees in `posterior_trees` represent the method's consensus concerning
phylogeny. This consensus will become clearer with future versions of our
method, when we will cluster these trees.


License
-------

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
