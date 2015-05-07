PhyloWGS input parser
=====================

Description
-----------
This parser can be used to create the `ssm_data.txt` and `cnv_data.txt` inputs
used by PhyloWGS.

The parser will create `ssm_data.txt` from a provided VCF file. Though VCF
itself is standardized, different variant callers differ in how they represent
read counts, and so minor adjustments must be made for each caller. Currently,
Sanger and MuTect VCFs are supported. To add support for your own flavour, you
may subclass `VariantParser` (see `SangerParser` and `MutectParser` for
examples).

For `cnv_data.txt`, the parser requires CNV calls. At the moment, only
[Battenberg](https://github.com/cancerit/cgpBattenberg) is supported. In the
future, we will likely support parsing
[TITAN](http://compbio.bccrc.ca/software/titan/) calls as well. Adding support
for a different CNV-caller is possible, but more difficult than supporting
other variant callers. Regardless, seeing how to generate the `d` value in
`cnv_data.txt`, which represents your confidence in the correctness of the
observed fraction of cells bearing the CNV, will be informative, as
establishing a proper value for this is difficult, given that it depends on the
size of the CNV, the read depth, the total copy-number change, and other
factors. For this code, please see the method `CnvFormatter._calc_total_reads`.

Installation
------------
The parser requires Python 2 and [PyVCF](https://pypi.python.org/pypi/PyVCF).
The latter can be installed via pip:

    pip2 install --user pyvcf

Usage
-----
    usage: create_phylowgs_inputs.py [-h] [-e ERROR_RATE] [-s SAMPLE_SIZE]
                                     [-b BATTENBERG] [--only-normal-cn]
                                     [--output-cnvs OUTPUT_CNVS]
                                     [--output-variants OUTPUT_VARIANTS]
                                     [-c CELLULARITY] -v {sanger,oncoscan,mutect}
                                     [--cnv-confidence CNV_CONFIDENCE]
                                     [--read-length READ_LENGTH] [--verbose]
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
      -b BATTENBERG, --battenberg BATTENBERG
                            Path to Battenberg-formatted list of CNVs (default:
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
                            somatic (default: 1.0)
      -v {sanger,oncoscan,mutect}, --variant-type {sanger,oncoscan,mutect}
                            Type of VCF file (default: None)
      --cnv-confidence CNV_CONFIDENCE
                            Confidence in CNVs. Set to < 1 to scale "d" values
                            used in CNV output file (default: 1.0)
      --read-length READ_LENGTH
                            Approximate length of reads. Used to calculate
                            confidence in CNV frequencies (default: 100)
      --verbose

Examples
--------
* Create ssm_data.txt from Sanger's PCAWG variant-calling
  pipeline, subsampling to 5000 mutations:

        ./create_phylowgs_inputs.py -s 5000 -v sanger sample.vcf

    * Note that no cnv_data.txt is created in this scenario, as no Battenberg
      data was provided. In this case, please provide an empty cnv_data.txt for
      PhyloWGS (e.g., created via `touch cnv_data.txt`).

* Create ssm_data.txt and cnv_data.txt from Sanger's PCAWG variant-calling
  pipeline, including all mutations:

        ./create_phylowgs_inputs.py -v sanger -b cnv_calls_from_battenberg.txt sample.vcf

* Only use variants in copy-number-normal regions:

        ./create_phylowgs_inputs.py -v sanger -b cnv_calls_from_battenberg.txt --only-normal-cn sample.vcf

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
  itself because of excessively high confidence in its population frequency.
