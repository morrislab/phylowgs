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

To see how to generate `ssm_data.txt` and `cnv_data.txt` from a VCF file and
[Battenberg](https://github.com/cancerit/cgpBattenberg) CNV file, please see
[the included parser](parser/).

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

`cnv_data.txt`: Note that if you are running without any CNVs, this file should
be empty. You can create the empty file via the command `touch cnv_data.txt`.

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

        g++ -o mh.o -O3 mh.cpp  util.cpp `gsl-config --cflags --libs`

3. Run PhyloWGS. Minimum invocation on sample data set:

        python2 evolve.py ssm_data.txt cnv_data.txt

  All options:

        usage: evolve.py [-h] [-b WRITE_BACKUPS_EVERY] [-k TOP_K_TREES]
                         [-f CLONAL_FREQS] [-s MCMC_SAMPLES] [-i MH_ITERATIONS]
                         [-r RANDOM_SEED]
                         ssm_file cnv_file

        Run PhyloWGS to infer subclonal composition from SSMs and CNVs

        positional arguments:
          ssm_file              File listing SSMs (simple somatic mutations, i.e.,
                                single nucleotide variants. For proper format, see
                                README.md.
          cnv_file              File listing CNVs (copy number variations). For proper
                                format, see README.md.

        optional arguments:
          -h, --help            show this help message and exit
          -b WRITE_BACKUPS_EVERY, --write-backups-every WRITE_BACKUPS_EVERY
                                Number of iterations to go between writing backups of
                                program state (default: 100)
          -k TOP_K_TREES, --top-k-trees TOP_K_TREES
                                Output file to save top-k trees in text format
                                (default: top_k_trees)
          -f CLONAL_FREQS, --clonal-freqs CLONAL_FREQS
                                Output file to save clonal frequencies (default:
                                clonalFrequencies)
          -s MCMC_SAMPLES, --mcmc-samples MCMC_SAMPLES
                                Number of MCMC samples (default: 2500)
          -i MH_ITERATIONS, --mh-iterations MH_ITERATIONS
                                Number of Metropolis-Hastings iterations (default:
                                5000)
          -r RANDOM_SEED, --random-seed RANDOM_SEED
                                Random seed for initializing MCMC sampler (default:
                                None)

4. Generate JSON results.

        mkdir test_results
        cd test_results
        # To work with viewer in Step 5, the naming conventions used here must be
        # followed.
        # "example_data" is simply the name by which you want your results to be identified.
        python2 /path/to/phylowgs/write_results.py example_data ../trees.zip example_data.summ.json.gz example_data.muts.json.gz example_data.mutass.zip
        cd ..

  All options:

        usage: write_results.py [-h] [--include-ssm-names] [--min-ssms MIN_SSMS]
                                dataset_name tree_file tree_summary_output
                                mutlist_output mutass_output

        Write JSON files describing trees

        positional arguments:
          dataset_name         Name identifying dataset
          tree_file            File containing sampled trees
          tree_summary_output  Output file for JSON-formatted tree summaries
          mutlist_output       Output file for JSON-formatted list of mutations
          mutass_output        Output file for JSON-formatted list of SSMs and CNVs
                               assigned to each subclone

        optional arguments:
          -h, --help           show this help message and exit
          --include-ssm-names  Include SSM names in output (which may be sensitive
                               data) (default: False)
          --min-ssms MIN_SSMS  Minimum number or percent of SSMs to retain a subclone
                               (default: 0.01)

5. View results.

        mv test_results /path/to/phylowgs/witness/data
        cd /path/to/phylowgs/witness
        gunzip data/*/*.gz
        python2 index_data.py
        python2 -m SimpleHTTPServer
        # Open http://127.0.0.1:8000 in your web browser. Note that, by
        # default, the server listens for connections from any host.


Resuming a previous PhyloWGS run
--------------------------------

If PhyloWGS is interrupted for whatever reason, you can resume your existing
run by simply running `evolve.py` from the same directory as the previous run,
without any command-line params:

    # Start initial run.
    python2 evolve.py ssm_data.txt cnv_data.txt

    # Hit CTRL+C to send SIGINT, halting run partway through.

    # Resume run:
    python2 evolve.py


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
