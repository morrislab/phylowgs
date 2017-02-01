#!/usr/bin/env python2
import argparse
from pwgsresults.result_munger import ResultMunger
from pwgsresults.result_loader import ResultLoader
from pwgsresults.json_writer import JsonWriter

def munge(treesummfn, mutlistfn, mutassfn):
  loader = ResultLoader(treesummfn, mutlistfn, mutassfn)
  dataset_name = loader.dataset_name
  treesumm, mutlist, mutass = loader.tree_summary, loader.mutlist, loader.load_all_mut_assignments_into_memory()

  munger = ResultMunger(treesumm, mutlist, mutass)
  munger.remove_superclones()
  munger.remove_polyclonal_trees()

  writer = JsonWriter(dataset_name)
  writer.write_summaries(treesumm, treesummfn)
  writer.write_mutlist(mutlist, mutlistfn)
  writer.write_mutass(mutass, mutassfn)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument('treesummfn',
    help='Output file for JSON-formatted tree summaries')
  parser.add_argument('mutlistfn',
    help='Output file for JSON-formatted list of mutations')
  parser.add_argument('mutassfn',
    help='Output file for JSON-formatted list of SSMs and CNVs assigned to each subclone')
  args = parser.parse_args()

  munge(args.treesummfn, args.mutlistfn, args.mutassfn)

if __name__ == '__main__':
	main()
