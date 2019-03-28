#!/usr/bin/env python2
import argparse
from pwgsresults.result_generator import ResultGenerator
from pwgsresults.result_munger import ResultMunger
from pwgsresults.json_writer import JsonWriter
from pwgsresults.tree_clusterer import TreeClusterer
from pwgsresults.ssm_analyser import SSM_Analyser
	

def main():
  parser = argparse.ArgumentParser(
    description='Write JSON files describing trees',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument('--include-ssm-names', dest='include_ssm_names', action='store_true',
    help='Include SSM names in output (which may be sensitive data)')
  parser.add_argument('--min-ssms', dest='min_ssms', type=float, default=0.01,
    help='Minimum number or percent of SSMs to retain a subclone')
  parser.add_argument('--clust-method', dest='clust_method', type=str, default="spectral",
    help='Which method to use when clustering sampled trees. Options are "gmm" and "spectral"(default)')
  parser.add_argument('dataset_name',
    help='Name identifying dataset')
  parser.add_argument('tree_file',
    help='File containing sampled trees')
  parser.add_argument('tree_summary_output',
    help='Output file for JSON-formatted tree summaries')
  parser.add_argument('mutlist_output',
    help='Output file for JSON-formatted list of mutations')
  parser.add_argument('mutass_output',
    help='Output file for JSON-formatted list of SSMs and CNVs assigned to each subclone')
  args = parser.parse_args()

  summaries, mutlist, mutass, params = ResultGenerator().generate(args.tree_file, args.include_ssm_names)

  munger = ResultMunger(summaries, mutlist, mutass)
  summaries, mutass = munger.remove_small_nodes(args.min_ssms)
  munger.remove_superclones()
  munger.remove_polyclonal_trees()
  
  clusters = TreeClusterer().find_clusters(summaries,args.clust_method) 
  SSM_Analyser().analyse(clusters, summaries, mutlist, mutass)

  writer = JsonWriter(args.dataset_name)
  writer.write_summaries(summaries, params, args.tree_summary_output, clusters)
  writer.write_mutlist(mutlist, args.mutlist_output)
  writer.write_mutass(mutass, args.mutass_output)

if __name__ == '__main__':
  main()
