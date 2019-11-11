cwlVersion: v1.0
class: CommandLineTool
label: Write Results
baseCommand: ["python", "/opt/phylowgs/write_results.py"]
requirements:
  - class: DockerRequirement
    dockerPull: smcheteval/phylowgs:0.1

inputs:
  dataset_name:
    type: string
    default: tumour
    inputBinding:
      position: 1

  tree_file:
    type: File
    inputBinding:
      position: 2

  tree_summary_output:
    type: string
    default: ./trees.json.gz
    inputBinding:
      position: 3

  mutlist_output:
    type: string
    default: ./mutations.json.gz
    inputBinding:
      position: 4

  mutass_output:
    type: string
    default: ./mutation_assignments.json.gz
    inputBinding:
      position: 5

outputs:
  summary_results:
    type: File
    outputBinding:
      glob: $(inputs.tree_summary_output)

  mutlist_results:
    type: File
    outputBinding:
      glob: $(inputs.mutlist_output)

  mutass_results:
    type: File
    outputBinding:
      glob: $(inputs.mutass_output)
