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
[TITAN](http://compbio.bccrc.ca/software/titan/) are supported. As the process of going from CNV calls to PhyloWGS input is complex -- the `a` and `d` values used by PhyloWGS for each CNV depend on the size of the CNV, the read depth, the total copy-number change, and other factors -- we perform CNV parsing as a two-step process:
  
    1. Run `parse_cnvs.py` to create the intermediate `cnvs.txt` file, which
       will contain for each CNV its chromosome, start and end coordinates, major
       and minor copy numbers, and clonal fraction (i.e., fraction of canceous
       cells containing the CNV, *not* the fraction of sample containing the CNV).

    2. Run `create_phylowgs_inputs.py` with the `--cnvs cnvs.txt` parameter.

Supporting other CNV callers is simply a matter of converting their results to
our intermediate `cnvs.txt` format. For an example of how this file should be
formatted, please see `cnvs.txt.example`.

Installation
------------
The parser requires Python 2 and [PyVCF](https://pypi.python.org/pypi/PyVCF).
The latter can be installed via pip:

    pip2 install --user pyvcf

Examples
--------
* Create ssm_data.txt from Sanger's PCAWG variant-calling
  pipeline, subsampling to 5000 mutations:

        ./create_phylowgs_inputs.py -s 5000 -v sanger sample.vcf

    * Note that an empty cnv_data.txt is created in this scenario, as no CNV
      data was provided.

* Create ssm_data.txt and cnv_data.txt from Sanger's PCAWG variant-calling
  pipeline, including all mutations, assuming 0.72 cellularity (i.e., sample purity):

        ./parse_cnvs.py -f battenberg cnv_calls_from_battenberg.txt 
        ./create_phylowgs_inputs.py -v sanger --cnvs cnvs.txt -c 0.72 sample.vcf

* Only use variants in copy-number-normal regions:

        ./parse_cnvs.py -f battenberg cnv_calls_from_battenberg.txt 
        ./create_phylowgs_inputs.py -v sanger -b cnv_calls_from_battenberg.txt -c 0.72 --only-normal-cn sample.vcf

Usage
-----
### CNV pre-parser

    usage: parse_cnvs.py [-h] -f {battenberg,titan}
                         [--cnv-output CNV_OUTPUT_FILENAME]
                         cnv_file

    Create CNV input file for parser from Battenberg or TITAN data

    positional arguments:
      cnv_file

    optional arguments:
      -h, --help            show this help message and exit
      -f {battenberg,titan}, --cnv-format {battenberg,titan}
                            Type of CNV input (default: None)
      --cnv-output CNV_OUTPUT_FILENAME
                            Output destination for parsed CNVs (default: cnvs.txt)

### Primary parser

    usage: create_phylowgs_inputs.py [-h] [-e ERROR_RATE] [-s SAMPLE_SIZE]
                                     [--cnvs CNV_FILE] [--only-normal-cn]
                                     [--output-cnvs OUTPUT_CNVS]
                                     [--output-variants OUTPUT_VARIANTS]
                                     [-c CELLULARITY] -v
                                     {sanger,mutect_pcawg,mutect_smchet,muse,dkfz,strelka,vardict}
                                     [--tumor-sample TUMOR_SAMPLE]
                                     [--cnv-confidence CNV_CONFIDENCE]
                                     [--read-length READ_LENGTH] [--verbose]
                                     [--muse-tier MUSE_TIER]
                                     vcf_file

    Create ssm_dat.txt and cnv_data.txt input files for PhyloWGS from VCF and CNV
    data.

    positional arguments:
      vcf_file

    optional arguments:
      -h, --help            show this help message and exit
      -e ERROR_RATE, --error-rate ERROR_RATE
                            Expected error rate of sequencing platform (default:
                            0.001)
      -s SAMPLE_SIZE, --sample-size SAMPLE_SIZE
                            Subsample SSMs to reduce PhyloWGS runtime (default:
                            None)
      --cnvs CNV_FILE       Path to CNV list created with parse_cnvs.py (default:
                            None)
      --only-normal-cn      Only output variants lying in normal CN regions. Do
                            not output CNV data directly. (default: False)
      --output-cnvs OUTPUT_CNVS
                            Output destination for CNVs (default: cnv_data.txt)
      --output-variants OUTPUT_VARIANTS
                            Output destination for variants (default:
                            ssm_data.txt)
      -c CELLULARITY, --cellularity CELLULARITY
                            Fraction of sample that is cancerous rather than
                            somatic. Used only for estimating CNV confidence -- if
                            no CNVs, need not specify argument. (default: 1.0)
      -v {sanger,mutect_pcawg,mutect_smchet,muse,dkfz,strelka,vardict}, --variant-type {sanger,mutect_pcawg,mutect_smchet,muse,dkfz,strelka,vardict}
                            Type of VCF file (default: None)
      --tumor-sample TUMOR_SAMPLE
                            Name of the tumor sample in the input VCF file.
                            Defaults to last sample if not specified. (default:
                            None)
      --cnv-confidence CNV_CONFIDENCE
                            Confidence in CNVs. Set to < 1 to scale "d" values
                            used in CNV output file (default: 1.0)
      --read-length READ_LENGTH
                            Approximate length of reads. Used to calculate
                            confidence in CNV frequencies (default: 100)
      --verbose
      --muse-tier MUSE_TIER
                            Maximum MuSE tier to include (default: 0)

Notes
-----
* Currently, limiting the number of variants used for phylogenetic
  reconstructions to 5000 is recommended to limit PhyloWGS' runtime. PhyloWGS
  runtime will scale linearly with the number of variants.

* Any variants on mitochondrial or sex chromosomes are ignored. Mitochondrial
  variants are weird, and sex-chromosome variants require knowing the patient's
  gender to generate proper allele frequencies.

* To permit CNVs to move more freely amongst populations in the sampled
  phylogenies, you can pass `--cnv-confidence [0.0, 1.0]`. This will scale the
  `d` column in `cnv_data.txt` by the given factor. A value of 1.0 leaves the
  `d` values unchanged; lower values will reduce every `d` value's magnitude by
  the specified factor, permitting CNVs to move more freely betwen populations,
  and reducing the probability that a CNV will be assigned a population unto
  itself because of excessively high confidence in its population frequency. To
  date, we've achieved our best results leaving this value at its default of
  1.0.
