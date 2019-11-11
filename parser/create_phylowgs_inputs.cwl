cwlVersion: v1.0
class: CommandLineTool
label: PhyloWGS Input Prep
baseCommand: ["python","/opt/phylowgs/parser/create_phylowgs_inputs.py"]
requirements:
  - class: DockerRequirement
    dockerPull: smcheteval/phylowgs:0.1

inputs:
  cnvs:
    type: File
    inputBinding:
      prefix: "--cnvs s1="
      separate: false

  output_cnvs:
    type: string
    default: cnv_data.txt
    inputBinding:
      prefix: --output-cnvs

  output_variants:
    type: string
    default: ssm_data.txt
    inputBinding:
      prefix: --output-variants

  vcf_type:
    type: string
    default: "s1=mutect_smchet"
    inputBinding:
      prefix: --vcf-type

  vcf_files:
    type: File
    inputBinding:
      prefix: "s1="
      separate: false
      position: 5

outputs:
  multievolve_cnvs:
    type: File
    outputBinding:
      glob: $(inputs.output_cnvs)

  multievolve_snvs:
    type: File
    outputBinding:
      glob: $(inputs.output_variants)
