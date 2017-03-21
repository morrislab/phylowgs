PhyloWGS input parser
=====================

Description
-----------
This parser can be used to create the `ssm_data.txt` and `cnv_data.txt` inputs
used by PhyloWGS.

The parser will create `ssm_data.txt` from a provided VCF file. Though VCF
itself is standardized, different variant callers differ in how they represent
read counts, and so minor adjustments must be made for each caller. A variety
of variant callers are already supported. To add support for your own flavour,
you may subclass `VariantParser` (see `SangerParser` or `MutectParser` for
examples).

For `cnv_data.txt`, the parser requires CNV calls.
[Battenberg](https://github.com/cancerit/cgpBattenberg) and
[TITAN](http://compbio.bccrc.ca/software/titan/) are supported. As the process
of going from CNV calls to PhyloWGS input is complex -- the `a` and `d` values
used by PhyloWGS for each CNV depend on the size of the CNV, the read depth,
the total copy-number change, and other factors -- we perform CNV parsing as a
two-step process:
  
  1. Run `parse_cnvs.py` to create the intermediate `cnvs.txt` file, which will
     contain for each CNV its chromosome, start and end coordinates, major and
     minor copy numbers, and cellular prevalence (i.e., fraction of cells in
     sample containing the CNV, *not* just the fraction of tumor cells containing
     the CNV).

  2. Run `create_phylowgs_inputs.py` with the `--cnvs <sample_name>=cnvs.txt` parameter.

`<sample_name>` is a string that uniquely identifies the given sample. This
isn't terribly useful when running on only a single sample; in such cases, you
can just name your sample `S1`, for example.
Supporting other CNV callers is simply a matter of converting their results to
our intermediate `cnvs.txt` format. For an example of how this file should be
formatted, please see `cnvs.txt.example`.

Please note that your listing of CNVs *must* list regions of normal copy number
as well. Any SSMs given in your VCF that fall outside regions listed in
`cnvs.txt` will be ignored by the parser, as it assumes your CNV caller could
not make a proper determination of copy number (whether normal or abnormal) for
unlisted regions, meaning that PhyloWGS will be unable to correct for
copy-number changes for the SSMs in question.

You can choose to run on three types of regions, specified via the `--regions`
parameter.

    * `--regions normal_and_abnormal_cn` (default): Include variants in clonal normal
      regions, clonal abnormal regions, and regions that have at most one
      subclonal abnormal state (e.g., a region that's a mix of (major, minor) =
      (2, 1) and (2, 2) will be rejected). Place CNAs on the phylogenetic tree.

    * `--regions normal_cn`: Include variants only in clonal normal regions
      (i.e., major and minor alleles are both 1). Do not place CNAs on the
      phylogenetic tree (i.e., an empty `cnv_data.txt` is output).

    * `--regions all`: *Recommended only if you have no data on CNAs, or if you
      are running PhyloWGS as part of the SMC-Het challenge*. Include every
      variant. This will include variants in regions with multiple subclonal
      abnormal states, as well as regions for which no copy number information
      is reported.

Installation
------------
The parser requires Python 2, NumPy, SciPy, and
[PyVCF](https://pypi.python.org/pypi/PyVCF). The last of these can be installed
via pip:

    pip2 install --user pyvcf

Examples
--------
Note that, in the following, the `sample1` name is arbitrary. As these are all
single-sample examples, the sample name isn't significant; nevertheless it must
be specified.

* Create ssm_data.txt from Sanger's PCAWG variant-calling
  pipeline, subsampling to 5000 mutations:

        ./create_phylowgs_inputs.py -s 5000 --vcf-type sample1=sanger sample1=sample.vcf

    * Note that an empty cnv_data.txt is created in this scenario, as no CNV
      data was provided.

* Create ssm_data.txt and cnv_data.txt from Sanger's PCAWG variant-calling
  pipeline, including all mutations, assuming 0.72 cellularity (i.e., sample purity):

        ./parse_cnvs.py -f battenberg -c 0.72 cnv_calls_from_battenberg.txt
        ./create_phylowgs_inputs.py --cnvs sample1=cnvs.txt --vcf-type sample1=sanger sample1=sample.vcf

* Only use variants in copy-number-normal regions:

        ./parse_cnvs.py -f battenberg -c 0.72 cnv_calls_from_battenberg.txt 
        ./create_phylowgs_inputs.py --cnvs  sample1=cnv_calls_from_battenberg.txt --only-normal-cn --vcf-type sample1=sanger sample1=sample.vcf

* Run with SSMs from VarDict and CNVs from TITAN, using cellularity of 0.81
  (which may be estimated from [1 - "Normal contamination estimate"] in TITAN's
  `*params.txt` output):

        ./parse_cnvs.py -f titan -c 0.81 cnv_calls_segs.txt
        ./create_phylowgs_inputs.py --cnvs sample1=cnvs.txt --vcf-type sample1=vardict sample1=sample.vcf

Multi-sample mode
-----------------
You can generate inputs from an arbitrary number of samples, each of which must have its own SSM and CNV calls. First, parse the CNV calls for each sample:

        ./parse_cnvs.py -f titan -c 0.81 --cnv-output cnvs1.txt samp1_segs.txt
        ./parse_cnvs.py -f titan -c 0.74 --cnv-output cnvs2.txt samp2_segs.txt

Now run the parser. Note that you must specify unique sample names for each
`--cnvs` option, each `--vcf-type` option, and each VCF file you pass. Here we
use sample names `S1` and `S2`, but they can be arbitrary strings. The order of
the VCF files given determines the order in which samples will appear during
post-processing; the order of the different `--cnvs` and `--vcf-type` options,
however, does not matter.

        ./create_phylowgs_inputs.py --cnvs S1=cnvs1.txt --cnvs S2=cnvs2.txt --vcf-type S1=vardict --vcf-type S2=vardict S1=variants1.vcf S2=variants2.vcf

Finally, run PhyloWGS. To retain in post-processing the sample names given to
the parser, include the `--params` option, passing the `params.json` file
generated by the parser.

        ../evolve.py --params params.json ssm_data.txt cnv_data.txt

Usage
-----
### CNV pre-parser

    usage: parse_cnvs.py [-h] -f {battenberg,titan} -c CELLULARITY
                         [--cnv-output CNV_OUTPUT_FILENAME]
                         cnv_file

    Create CNV input file for parser from Battenberg or TITAN data

    positional arguments:
      cnv_file

    optional arguments:
      -h, --help            show this help message and exit
      -f {battenberg,titan}, --cnv-format {battenberg,titan}
                            Type of CNV input (default: None)
      -c CELLULARITY, --cellularity CELLULARITY
                            Fraction of sample that is cancerous rather than
                            somatic. Used only for estimating CNV confidence -- if
                            no CNVs, need not specify argument. (default: None)
      --cnv-output CNV_OUTPUT_FILENAME
                            Output destination for parsed CNVs (default: cnvs.txt)
### Primary parser

    usage: create_phylowgs_inputs.py [-h] --vcf-type VCF_TYPES [-e ERROR_RATE]
                                     [--missing-variant-confidence MISSING_VARIANT_CONFIDENCE]
                                     [-s SAMPLE_SIZE] [-P PRIORITY_SSM_FILENAME]
                                     [--cnvs CNV_FILES]
                                     [--regions {normal_cn,normal_and_abnormal_cn,all}]
                                     [--output-cnvs OUTPUT_CNVS]
                                     [--output-variants OUTPUT_VARIANTS]
                                     [--output-params OUTPUT_PARAMS]
                                     [--tumor-sample TUMOR_SAMPLE]
                                     [--cnv-confidence CNV_CONFIDENCE]
                                     [--read-length READ_LENGTH]
                                     [--muse-tier MUSE_TIER]
                                     [--nonsubsampled-variants OUTPUT_NONSUBSAMPLED_VARIANTS]
                                     [--nonsubsampled-variants-cnvs OUTPUT_NONSUBSAMPLED_VARIANTS_CNVS]
                                     [--sex {auto,male,female}] [--verbose]
                                     vcf_files [vcf_files ...]

    Create ssm_dat.txt and cnv_data.txt input files for PhyloWGS from VCF and CNV
    data.

    positional arguments:
      vcf_files             Path to VCF file for each sample. Specified as
                            <sample>=<VCF path>.

    optional arguments:
      -h, --help            show this help message and exit
      --vcf-type VCF_TYPES  Type of VCF file for each sample, specified as
                            <sample>=<vcf_type>. Valid VCF types are strelka,mutec
                            t_pcawg,dkfz,muse,vardict,mutect_smchet,mutect_tcga,sa
                            nger,pcawg_consensus. (default: None)
      -e ERROR_RATE, --error-rate ERROR_RATE
                            Expected error rate of sequencing platform (default:
                            0.001)
      --missing-variant-confidence MISSING_VARIANT_CONFIDENCE
                            Confidence in range [0, 1] that SSMs missing from a
                            sample are indeed not present in that sample (default:
                            1.0)
      -s SAMPLE_SIZE, --sample-size SAMPLE_SIZE
                            Subsample SSMs to reduce PhyloWGS runtime (default:
                            None)
      -P PRIORITY_SSM_FILENAME, --priority-ssms PRIORITY_SSM_FILENAME
                            File containing newline-separated list of SSMs in
                            "<chr>_<locus>" format to prioritize for inclusion
                            (default: None)
      --cnvs CNV_FILES      Path to CNV file created with parse_cnvs.py for each
                            sample. Specified as <sample>=<CNV path>. (default:
                            None)
      --regions {normal_cn,normal_and_abnormal_cn,all}
                            Which regions to use variants from. Refer to the
                            parser README for more details. (default:
                            normal_and_abnormal_cn)
      --output-cnvs OUTPUT_CNVS
                            Output destination for CNVs (default: cnv_data.txt)
      --output-variants OUTPUT_VARIANTS
                            Output destination for variants (default:
                            ssm_data.txt)
      --output-params OUTPUT_PARAMS
                            Output destination for run parameters (default:
                            params.json)
      --tumor-sample TUMOR_SAMPLE
                            Name of the tumor sample in the input VCF file.
                            Defaults to last sample if not specified. (default:
                            None)
      --cnv-confidence CNV_CONFIDENCE
                            Confidence in CNVs. Set to < 1 to scale "d" values
                            used in CNV output file (default: 0.5)
      --read-length READ_LENGTH
                            Approximate length of reads. Used to calculate
                            confidence in CNV frequencies (default: 100)
      --muse-tier MUSE_TIER
                            Maximum MuSE tier to include (default: 0)
      --nonsubsampled-variants OUTPUT_NONSUBSAMPLED_VARIANTS
                            If subsampling, write nonsubsampled variants to
                            separate file, in addition to subsampled variants
                            (default: None)
      --nonsubsampled-variants-cnvs OUTPUT_NONSUBSAMPLED_VARIANTS_CNVS
                            If subsampling, write CNVs for nonsubsampled variants
                            to separate file (default: None)
      --sex {auto,male,female}
                            Sex of patient. Used to adjust expected variant
                            frequencies on sex chromosomes. If auto, patient is
                            set to male if any variants are provided on the Y
                            chromosome, and female otherwise. (default: auto)
      --verbose

Notes
-----
* Currently, limiting the number of variants used for phylogenetic
  reconstructions to 5000 is recommended to limit PhyloWGS' runtime. PhyloWGS
  runtime will scale linearly with the number of variants.

* Any variants on the mitochondrial chromosome are ignored. Mitochondrial
  variants are weird.

* To permit CNVs to move more freely amongst populations in the sampled
  phylogenies, you can pass `--cnv-confidence [0.0, 1.0]`. This will scale the
  `d` column in `cnv_data.txt` by the given factor. A value of 1.0 leaves the
  `d` values unchanged; lower values will reduce every `d` value's magnitude by
  the specified factor, permitting CNVs to move more freely betwen populations,
  and reducing the probability that a CNV will be assigned a population unto
  itself because of excessively high confidence in its population frequency. To
  date, we've achieved our best results leaving this value at its default of
  0.5.
