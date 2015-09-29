#!/usr/bin/env python2
import argparse
import json
import gzip
import zipfile
from pwgsresults.result_generator import ResultGenerator
from pwgsresults.result_munger import ResultMunger

class JsonWriter(object):
  def __init__(self, dataset_name):
    self._dataset_name = dataset_name

  def write_mutlist(self, mutlist, mutlist_outfn):
    with gzip.GzipFile(mutlist_outfn, 'w') as mutf:
      mutlist['dataset_name'] = self._dataset_name
      json.dump(mutlist, mutf)

  def write_summaries(self, summaries, summaries_outfn):
    to_dump = {
      'dataset_name': self._dataset_name,
      'trees': summaries,
    }
    with gzip.GzipFile(summaries_outfn, 'w') as summf:
      json.dump(to_dump, summf)

  def write_mutass(self, mutass, mutass_outfn):
    with zipfile.ZipFile(mutass_outfn, 'w', compression=zipfile.ZIP_DEFLATED) as muts_file:
      for tree_idx, tree_mutass in mutass.items():
        to_dump = {
          'mut_assignments': tree_mutass,
          'dataset_name': self._dataset_name
        }
        muts_file.writestr('%s.json' % tree_idx, json.dumps(to_dump))

def main():
  parser = argparse.ArgumentParser(
    description='Write JSON files describing trees',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument('--include-ssm-names', dest='include_ssm_names', action='store_true',
    help='Include SSM names in output (which may be sensitive data)')
  parser.add_argument('--min-ssms', dest='min_ssms', type=float, default=0.01,
    help='Minimum number or percent of SSMs to retain a subclone')
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

  summaries, mutlist, mutass = ResultGenerator().generate(args.tree_file, args.include_ssm_names)

  munger = ResultMunger(summaries, mutlist, mutass, args.min_ssms)
  summaries, mutass = munger.munge()

  writer = JsonWriter(args.dataset_name)
  writer.write_summaries(summaries, args.tree_summary_output)
  writer.write_mutlist(mutlist, args.mutlist_output)
  writer.write_mutass(mutass, args.mutass_output)

if __name__ == '__main__':
  main()
