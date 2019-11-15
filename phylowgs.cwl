cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement

inputs:
  ssmFile:
    type: File
  cnvFile:
    type: File
  cellFxn:
    type: float

steps:
  parse:
    run: ./parser/parse_cnvs.cwl
    in:
      cellularity: cellFxn
      cnv_file: cnvFile
    out:
      - parser_output

  prep_inputs:
    run: ./parser/create_phylowgs_inputs.cwl
    in:
      cnvs: parse/parser_output
      vcf_files: ssmFile
    out:
      - multievolve_snvs
      - multievolve_cnvs

  multievolve:
    run: ./multievolve.cwl
    in:
      ssms: prep_inputs/multievolve_snvs
      cnvs: prep_inputs/multievolve_cnvs
    out:
      - tree_file

  write_results:
    run: ./write_results.cwl
    in: 
      tree_file: multievolve/tree_file
    out:
      - summary_results
      - mutlist_results
      - mutass_results

  write_report:
    run: ../smchet-challenge/create-smchet-report/write_report.cwl
    in:
      tree_summary: write_results/summary_results
      mutation_list: write_results/mutlist_results
      mutation_assignment: write_results/mutass_results
    out:
      - cellularity
      - population
      - proportion
      - cluster_assignment
      - cocluster_assignment

outputs:
  cellularity_predfile:
    type: File
    outputSource: write_report/cellularity
  population_predfile:
    type: File
    outputSource: write_report/population
  proportion_predfile:
    type: File
    outputSource: write_report/proportion
  cluster_assignment_predfile:
    type: File
    outputSource: write_report/cluster_assignment
  cocluster_assignment_predfile:
    type: File
    outputSource: write_report/cocluster_assignment
