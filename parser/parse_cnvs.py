from __future__ import print_function
import argparse
import csv
import sys
from collections import defaultdict

def chrom_key(chrom):
  chrom = chrom.lower()
  if chrom == 'x':
    chrom = 100
  elif chrom == 'y':
    chrom = 101
  elif chrom.isdigit():
    chrom = int(chrom)
  else:
    chrom = 999
  return chrom

class CopyNumberWriter(object):
  def __init__(self, cn_output_fn):
    self._cn_output_fn = cn_output_fn
    # cellular_prevalence represents fraction of *all* cells that are affected
    # by CNVs, *not* just tumor cells.
    self._keys = ('chromosome', 'start', 'end', 'copy_number', 'minor_cn', 'major_cn', 'cellular_prevalence')

  def _write_header(self):
    self._cn_output.write('\t'.join(self._keys) + '\n')

  def _write_cn_record(self, region):
    vals = [str(region[k]) for k in self._keys]
    self._cn_output.write('\t'.join(vals) + '\n')

  def write_cnvs(self, cn_regions):
    self._cn_output = open(self._cn_output_fn, 'w')
    self._write_header()

    chroms = sorted(cn_regions.keys(), key = chrom_key)
    for chrom in chroms:
      chrom_regions = cn_regions[chrom]
      chrom_regions.sort(key = lambda r: r['start'])
      for region in chrom_regions:
        # Insert chromosome into record, as including it originally would have
        # duplicated the dictionary key corresponding to per-chromosome CNVs.
        region['chromosome'] = chrom
        region['copy_number'] = region['major_cn'] + region['minor_cn']
        self._write_cn_record(region)

    self._cn_output.close()

class CnvParser(object):
  def parse(self):
    raise Exception('Not implemented')

class TitanParser(CnvParser):
  def __init__(self, titan_filename, cellularity):
    self._titan_filename = titan_filename
    self._cellularity = cellularity

  def parse(self):
    cn_regions = defaultdict(list)

    with open(self._titan_filename) as titanf:
      reader = csv.DictReader(titanf, delimiter='\t')
      for record in reader:
        chrom = record['Chromosome'].lower()
        cnv = {}
        cnv['start'] = int(record['Start_Position(bp)'])
        cnv['end'] = int(record['End_Position(bp)'])
        cnv['major_cn'] = int(record['MajorCN'])
        cnv['minor_cn'] = int(record['MinorCN'])

        clonal_freq = record['Clonal_Frequency']
        if clonal_freq == 'NA':
          cnv['cellular_prevalence'] = self._cellularity
        else:
          cnv['cellular_prevalence'] = float(clonal_freq) * self._cellularity

        cn_regions[chrom].append(cnv)

    return cn_regions

class BattenbergParser(CnvParser):
  def __init__(self, bb_filename, cellularity):
    self._bb_filename = bb_filename
    self._cellularity = cellularity
    # Used by SMC-Het parser, which has fields shifted by 1.
    self._field_offset = 0

  def _compute_cn(self, cnv1, cnv2):
    '''
    This code isn't used, but is retained for reference.
    '''
    cn1 = (cnv1['nmaj'] + cnv1['nmin']) * cnv1['frac']
    if cnv2:
      cn2 = (cnv2['nmaj'] + cnv2['nmin']) * cnv2['frac']
    else:
      cn2 = 0
    total_cn = cn1 + cn2
    return total_cn

  def parse(self):
    cn_regions = defaultdict(list)
    pval_threshold = 0.05

    with open(self._bb_filename) as bbf:
      header = bbf.next()
      for line in bbf:
        fields = line.strip().split()
        chrom = fields[1 + self._field_offset].lower()
        start = int(fields[2 + self._field_offset])
        end = int(fields[3 + self._field_offset])
        pval = float(fields[5 + self._field_offset])

        cnv1 = {}
        cnv1['start'] = start
        cnv1['end'] = end
        cnv1['major_cn'] = int(fields[8 + self._field_offset])
        cnv1['minor_cn'] = int(fields[9 + self._field_offset])
        cnv1['cellular_prevalence'] = float(fields[10 + self._field_offset]) * self._cellularity

        cnv2 = None
        # Stefan's comment on p values: The p-values correspond "to whether a
        # segment should be clonal or subclonal copynumber. We first fit a
        # clonal copynumber profile for the whole sample and then perform a
        # simple two-sided t-test twhere the null hypothesis is: A particular
        # segment is clonal. And the alternative: It is subclonal."
        #
        # Thus: if t-test falls below significance threshold, we push cnv1 to
        # clonal frequency.
        if pval <= pval_threshold:
          cnv2 = {}
          cnv2['start'] = start
          cnv2['end'] = end
          cnv2['major_cn'] = int(fields[11 + self._field_offset])
          cnv2['minor_cn'] = int(fields[12 + self._field_offset])
          cnv2['cellular_prevalence'] = float(fields[13 + self._field_offset]) * self._cellularity
        else:
          cnv1['cellular_prevalence'] = self._cellularity

        cn_regions[chrom].append(cnv1)
        if cnv2 is not None:
          cn_regions[chrom].append(cnv2)
    return cn_regions

class BattenbergSmchetParser(BattenbergParser):
  def __init__(self, bb_filename, cellularity):
    super(BattenbergSmchetParser, self).__init__(bb_filename, cellularity)
    # SMC-Het Battenberg files lack the initial index column.
    self._field_offset = -1

def restricted_float(x):
  x = float(x)
  if x < 0.0 or x > 1.0:
    raise argparse.ArgumentTypeError('%r not in range [0.0, 1.0]' % x)
  return x

def main():
  parser = argparse.ArgumentParser(
    description='Create CNV input file for parser from Battenberg or TITAN data',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('-f', '--cnv-format', dest='input_type', required=True, choices=('battenberg', 'battenberg-smchet', 'titan'),
    help='Type of CNV input')
  parser.add_argument('-c', '--cellularity', dest='cellularity', type=float, required=True,
    help='Fraction of sample that is cancerous rather than somatic. Used only for estimating CNV confidence -- if no CNVs, need not specify argument.')
  parser.add_argument('--cnv-output', dest='cnv_output_filename', default='cnvs.txt',
    help='Output destination for parsed CNVs')
  parser.add_argument('cnv_file')
  args = parser.parse_args()

  if args.cellularity > 1.0:
    print('Cellularity for %s is %s. Setting to 1.0.' % (args.cnv_file, args.cellularity), file=sys.stderr)
    cellularity = 1.0
  else:
    cellularity = args.cellularity


  if args.input_type == 'battenberg':
    parser = BattenbergParser(args.cnv_file, cellularity)
  elif args.input_type == 'battenberg-smchet':
    parser = BattenbergSmchetParser(args.cnv_file, cellularity)
  elif args.input_type == 'titan':
    parser = TitanParser(args.cnv_file, cellularity)
  else:
    raise Exception('Unknown input type')

  writer = CopyNumberWriter(args.cnv_output_filename)
  regions = parser.parse()
  writer.write_cnvs(regions)

if __name__ == '__main__':
  main()
