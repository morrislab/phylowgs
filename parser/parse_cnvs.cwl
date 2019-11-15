cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["python","/opt/phylowgs/parser/parse_cnvs.py"]
requirements:
  - class: DockerRequirement
    dockerPull: smcheteval/phylowgs:0.1

inputs:
  cnv_format:
    type: string
    default: battenberg-smchet
    inputBinding:
      prefix: --cnv-format

  cellularity:
    type: float
    default: None
    inputBinding:
      prefix: --cellularity

  cnv_output:
    type: string
    default: cnvs.txt
    inputBinding:
      prefix: --cnv-output

  cnv_file:
    type: File
    inputBinding:
      position: 4

outputs:
  parser_output:
    type: File
    outputBinding:
      glob: $(inputs.cnv_output)
