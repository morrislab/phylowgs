cwlVersion: v1.0
class: CommandLineTool
label: Multievolve
baseCommand: ["python", "/opt/phylowgs/multievolve.py"]
requirements:
  - class: DockerRequirement
    dockerPull: smcheteval/phylowgs:0.1

inputs:
  num_chains:
    type: int
    default: 16
    inputBinding:
      prefix: --num-chains

  random_seeds:
    type: int[]
    default: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    inputBinding:
      prefix: --random-seeds

  chain_inclusion_factor:
    type: float
    default: 1.1
    inputBinding:
      prefix: --chain-inclusion-factor

  mcmc:
    type: int
    default: 5000
    inputBinding:
      prefix: --mcmc-samples

  burnin:
    type: int
    default: 2000
    inputBinding:
      prefix: --burnin-samples

  ssms:
    type: File
    inputBinding:
      prefix: --ssms

  cnvs:
    type: File
    inputBinding:
      prefix: --cnvs

  output_dir:
    type: string
    default: ./
    inputBinding:
      prefix: --output-dir

outputs:
  tree_file:
    type: File
    outputBinding:
      glob: trees.zip

