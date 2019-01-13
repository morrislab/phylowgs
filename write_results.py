#!/usr/bin/env python2
import argparse
from pwgsresults.result_generator import ResultGenerator
from pwgsresults.result_munger import ResultMunger
from pwgsresults.json_writer import JsonWriter

def restricted_float(X):
  X = float(X)
  L, U = 0, 1
  if not (L <= X <= U):
    raise argparse.ArgumentTypeError('%s outside range [%s, %s]' % (X, L, U))
  return X

def main():
  parser = argparse.ArgumentParser(
    description='Write JSON files describing trees',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument('--include-ssm-names', dest='include_ssm_names', action='store_true',
    help='Include SSM names in output (which may be sensitive data)')
  parser.add_argument('--min-ssms', dest='min_ssms', type=float, default=0.01,
    help='Minimum number or percent of SSMs to retain a subclone')
  parser.add_argument('--allow-multiprimary', dest='allow_multiprimary', action='store_true',
    help='Minimum number or percent of SSMs to retain a subclone')
  parser.add_argument('--include-multiprimary', dest='include_multiprimary', action='store_true',
    help='Whether to include multiprimary trees in result')
  parser.add_argument('--max-multiprimary', dest='max_multiprimary', type=restricted_float, default=0.8,
    help='Maximum proportion of trees that may be multiprimary if ' \
    '--include=multiprimary=False. In that case, An exception will be thrown if ' \
    'the proportion of multiprimary trees exceeds this value.')
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
  if not args.include_multiprimary:
    munger.remove_multiprimary_trees(args.max_multiprimary)

  writer = JsonWriter(args.dataset_name)
  writer.write_summaries(summaries, params, args.tree_summary_output)
  writer.write_mutlist(mutlist, args.mutlist_output)
  writer.write_mutass(mutass, args.mutass_output)

if __name__ == '__main__':
  main()
