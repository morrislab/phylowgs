#!/usr/bin/env python2
from __future__ import print_function
import argparse

def split(ssm_file):
  with open(ssm_file) as ssmf:
    header = next(ssmf).strip()
    for line in ssmf:
      line = line.strip()
      fields = line.split('\t')
      ssm_id = fields[0]
      with open('%s.ssm' % ssm_id, 'w') as outf:
        print(header, file=outf)
        print(line, file=outf)

def main():
  parser = argparse.ArgumentParser(
    description='Create separate SSM file for each SSM listed in a single SSM file',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('ssm_file', help='SSM file to split')

  args = parser.parse_args()
  split(args.ssm_file)

if __name__ == '__main__':
  main()
