PhyloWGS input parser
=====================

Description
-----------
This parser can be used to create the `ssm_data.txt` and `cnv_data.txt` inputs used by PhyloWGS.

Usage
-----
    usage: create_phylowgs_inputs.py [-h] [-e ERROR_RATE] [-s SAMPLE_SIZE]
                                     [-b BATTENBERG] [--only-normal-cn]
                                     [--output-cnvs OUTPUT_CNVS]
                                     [--output-variants OUTPUT_VARIANTS]
                                     [-c CELLULARITY] -v {sanger,oncoscan,mutect}
                                     [--verbose]
                                     vcf_file

    Create SSM input file for PhyloWGS from VCF and CNV data

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
