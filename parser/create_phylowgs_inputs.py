#!/usr/bin/env python2
from __future__ import print_function

# Requires PyVCF. To install: pip2 install pyvcf
import vcf
import argparse
import csv
from collections import defaultdict, namedtuple, OrderedDict
import random
import sys
import numpy as np
import numpy.ma as ma
from scipy.stats.mstats import gmean

class ReadCountsUnavailableError(Exception):
  pass

class VariantParser(object):
  def __init__(self):
    # Child classes must give the following variables sensible values in
    # constructor so that list_variants() works subsequently.
    self._cnvs = None
    self._vcf_filename = None

  def list_variants(self):
    variants = self._filter(self._vcf_filename)
    variants_and_reads = []
    for variant in variants:
      try:
        ref_reads, total_reads = self._calc_read_counts(variant)
      except ReadCountsUnavailableError as exc:
        continue
      variants_and_reads.append((variant, ref_reads, total_reads))
    return variants_and_reads

  def _calc_read_counts(self, variant):
    raise Exception('Not implemented -- use child class')

  def _parse_vcf(self, vcf_filename):
    vcfr = vcf.Reader(filename=vcf_filename)
    records = []
    for variant in vcfr:
      variant.CHROM = variant.CHROM.upper()
      # Some VCF dialects prepend "chr", some don't. Remove the prefix to
      # standardize.
      if variant.CHROM.startswith('CHR'):
        variant.CHROM = variant.CHROM[3:]
      records.append(variant)
    return records

  def _does_variant_pass_filters(self, variant):
    if variant.FILTER is None:
      return True
    if len(variant.FILTER) > 0:
      # Variant failed one or more filters.
      return False
    return True

  def _filter(self, vcf_filename):
    variants = []

    all_variants = self._parse_vcf(vcf_filename)

    for variant in all_variants:
      if not is_good_chrom(variant.CHROM):
        continue
      if not self._does_variant_pass_filters(variant):
        continue
      variants.append(variant)
    return variants

  def _get_tumor_index(self, variant, tumor_sample=None):
    """Find the index of the tumor sample.

    Currently hardcodes tumour sample as the last column if name not specified.
    Might not always be true
    """
    if self._tumor_sample:
      tumor_is = [i for i, s in enumerate(variant.samples) if s.sample == tumor_sample]
      assert len(tumor_is) == 1, "Did not find tumor name %s in samples" % tumor_sample
      return tumor_is[0]
    else:
      # Don't make this -1, as some code assumes it will be >= 0.
      return len(variant.samples) - 1

class SangerParser(VariantParser):
  '''
  Works with PCAWG variant calls from the Sanger.
  '''
  def __init__(self, vcf_filename, tumor_sample=None):
    self._vcf_filename = vcf_filename
    self._tumor_sample = tumor_sample

  def _find_ref_and_variant_nt(self, variant):
    assert len(variant.REF) == len(variant.ALT) == 1
    return (str(variant.REF[0]), str(variant.ALT[0]))

  def _calc_read_counts(self, variant):
    normal = variant.genotype('NORMAL')
    tumor = variant.genotype('TUMOUR')

    reference_nt, variant_nt = self._find_ref_and_variant_nt(variant)
    tumor_reads = {
      'forward': {
        'A': int(tumor['FAZ']),
        'C': int(tumor['FCZ']),
        'G': int(tumor['FGZ']),
        'T': int(tumor['FTZ']),
      },
      'reverse': {
        'A': int(tumor['RAZ']),
        'C': int(tumor['RCZ']),
        'G': int(tumor['RGZ']),
        'T': int(tumor['RTZ']),
      },
    }

    ref_reads = tumor_reads['forward'][reference_nt] + tumor_reads['reverse'][reference_nt]
    # For now, variant reads are defined as only the non-reference nucleotide in
    # the inferred tumor SNP. We ignore reads of a third or fourth base.
    variant_reads = tumor_reads['forward'][variant_nt] + tumor_reads['reverse'][variant_nt]
    total_reads = ref_reads + variant_reads

    return (ref_reads, total_reads)

class PcawgConsensusParser(VariantParser):
  def __init__(self, vcf_filename, tumor_sample=None):
    self._vcf_filename = vcf_filename
    self._tumor_sample = tumor_sample

  def _find_ref_and_variant_nt(self, variant):
    assert len(variant.REF) == len(variant.ALT) == 1
    return (str(variant.REF[0]), str(variant.ALT[0]))

  def _calc_read_counts(self, variant):
    if not ('t_alt_count' in variant.INFO and 't_ref_count' in variant.INFO):
      raise ReadCountsUnavailableError()
    assert len(variant.INFO['t_alt_count']) == len(variant.INFO['t_ref_count']) == 1

    alt_reads = int(variant.INFO['t_alt_count'][0])
    ref_reads = int(variant.INFO['t_ref_count'][0])
    total_reads = alt_reads + ref_reads
    # Some variants havezero alt and ref reads.
    if total_reads == 0:
      raise ReadCountsUnavailableError()
    return (ref_reads, total_reads)

class MuseParser(VariantParser):
  def __init__(self, vcf_filename, tier=0, tumor_sample=None):
    self._vcf_filename = vcf_filename
    self._tier = tier
    self._tumor_sample = tumor_sample

  def _get_normal_genotype(self, variant):
    tumor_i = self._get_tumor_index(variant, self._tumor_sample)
    assert tumor_i in (0, 1), 'Tumor index %s is not 0 or 1' % tumor_i
    normal_i = 1 - tumor_i
    return set([int(t) for t in variant.samples[normal_i]['GT'].split('/')])

  def _calc_read_counts(self, variant):
    normal_gt = self._get_normal_genotype(variant)
    assert len(normal_gt) == 1
    normal_gt = normal_gt.pop()

    tumor_i = self._get_tumor_index(variant, self._tumor_sample)
    total_reads = int(variant.samples[tumor_i]['DP'])
    ref_reads = int(variant.samples[tumor_i]['AD'][normal_gt])

    return (ref_reads, total_reads)

  def _does_variant_pass_filters(self, variant):
    # Ignore heterozygous normal variants.
    if len(self._get_normal_genotype(variant)) != 1:
      return False
    if variant.FILTER is None or len(variant.FILTER) == 0:
      return True
    if int(variant.FILTER[0][-1]) <= self._tier:
      # Variant failed one or more filters, but we still accept it.
      return True
    return False
    
class StrelkaParser(VariantParser):
  def __init__(self, vcf_filename, tumor_sample=None):
    self._vcf_filename = vcf_filename
    self._tumor_sample = tumor_sample    

  def _does_variant_pass_filters(self, variant):
    # Strelka outputs two files one for SNPs, the other for InDels
    # For now only deal with SNP file from Strelka
    if variant.is_snp:
      if variant.FILTER is None or len(variant.FILTER) == 0: 
        return True
    return False

  def _calc_read_counts(self, variant):
    alt = variant.ALT[0]
    tumor_i = self._get_tumor_index(variant, self._tumor_sample)
    total_reads = int(variant.samples[tumor_i]['DP'])

    if alt is None:
      total_reads = 0
      variant_reads = 0
    else:
      variant_reads = int(getattr(variant.samples[tumor_i].data, str(alt)+'U')[0])

    ref_reads = total_reads - variant_reads
    return (ref_reads, total_reads)

class MutectTcgaParser(VariantParser):
  def __init__(self, vcf_filename, tumor_sample=None):
    self._vcf_filename = vcf_filename
    self._tumor_sample = tumor_sample

  def _calc_read_counts(self, variant):
    tumor_i = self._get_tumor_index(variant, self._tumor_sample)
    # TD: Tumor allelic depths for the ref and alt alleles in the order listed
    ref_reads, variant_reads = variant.samples[tumor_i]['TD']
    total_reads = ref_reads + variant_reads
    return (ref_reads, total_reads)

class MutectPcawgParser(VariantParser):
  def __init__(self, vcf_filename, tumor_sample=None):
    self._vcf_filename = vcf_filename
    self._tumor_sample = tumor_sample

  def _calc_read_counts(self, variant):
    tumor_i = self._get_tumor_index(variant, self._tumor_sample)
    ref_reads = int(variant.samples[tumor_i].data.ref_count)
    variant_reads = int(variant.samples[tumor_i].data.alt_count)
    total_reads = ref_reads + variant_reads

    return (ref_reads, total_reads)

class MutectSmchetParser(VariantParser):
  def __init__(self, vcf_filename, tumor_sample=None):
    self._vcf_filename = vcf_filename
    self._tumor_sample = tumor_sample

  def _calc_read_counts(self, variant):
    tumor_i = self._get_tumor_index(variant, self._tumor_sample)
    ref_reads = int(variant.samples[tumor_i]['AD'][0])
    variant_reads = int(variant.samples[tumor_i]['AD'][1])
    total_reads = ref_reads + variant_reads

    return (ref_reads, total_reads)

class VarDictParser(MutectSmchetParser):
  """Support VarDict somatic variant caller.

  https://github.com/AstraZeneca-NGS/VarDictJava
  https://github.com/AstraZeneca-NGS/VarDict

  Uses the same read-extraction logic as MuTect (SMC-Het).
  """
  pass

class DKFZParser(VariantParser):
  def __init__(self, vcf_filename, tumor_sample=None):
    self._vcf_filename = vcf_filename
    self._tumor_sample = tumor_sample

  def _calc_read_counts(self, variant):
    # This doesn't handle multisample correctly, as I don't know how to get the
    # DP4 attribute on multiple DKFZ samples currently.
    for_ref_reads = int(variant.INFO['DP4'][0])
    back_ref_reads = int(variant.INFO['DP4'][1])
    for_variant_reads = int(variant.INFO['DP4'][2])
    back_variant_reads = int(variant.INFO['DP4'][3])
    ref_reads = for_ref_reads + back_ref_reads
    var_reads = for_variant_reads + back_variant_reads
    total_reads = ref_reads + var_reads

    return (ref_reads, total_reads)

class CnvFormatter(object):
  def __init__(self, cnv_confidence, read_depth, read_length, sampidxs):
    self._cnv_confidence = cnv_confidence
    self._read_depth = read_depth
    self._read_length = read_length
    self._sampidxs = sampidxs

  def _max_reads(self, sampidx):
    return 1e6 * self._read_depth[sampidx]

  def _find_overlapping_variants(self, chrom, cnv, variants):
    overlapping = []

    start = cnv['start']
    end = cnv['end']
    for variant in variants:
      if chrom.upper() == variant['chrom'].upper():
        if start <= variant['pos'] <= end:
          overlapping.append(variant['ssm_id'])
    return overlapping

  def _calc_ref_reads(self, cellular_prev, total_reads):
    ref_reads = np.zeros(len(self._sampidxs))
    for sampidx in self._sampidxs:
      vaf = cellular_prev[sampidx] / 2
      ref_reads[sampidx] = int((1 - vaf) * total_reads[sampidx])
    return ref_reads

  def _calc_total_reads(self, cellular_prev, locus_start, locus_end, major_cn, minor_cn):
    total_reads = np.zeros(len(self._sampidxs))
    assert len(set(major_cn)) == len(set(minor_cn)) == 1
    total_cn = major_cn[0] + minor_cn[0]

    for sampidx in self._sampidxs:
      # Proportion of all cells carrying CNV.
      P = cellular_prev[sampidx]
      if total_cn == 2:
        # If no net change in copy number -- e.g., because (major, minor) went
        # from (1, 1) to (2, 0) -- force the delta_cn to be 1.
        delta_cn = 1.
        no_net_change = True
      else:
        delta_cn = float(total_cn - 2)
        no_net_change = False

      region_length = locus_end - locus_start + 1
      fn = (self._read_depth[sampidx] * region_length) / self._read_length[sampidx]

      # This is a hack to prevent division by zero (when delta_cn = -2). Its
      # effect will be to make d large.
      if P > 0.98:
        P = 0.98
      if P < 0.02:
        P = 0.02

      d = (delta_cn**2 / 4) * (fn * P * (2 - P)) / (1 + (abs(delta_cn)  * P) / 2)

      if no_net_change:
        # If no net change in CN occurred, the estimate was just based on BAFs,
        # meaning we have lower confidence in it. Indicate this lack of
        # confidence via d by multiplying it by (read length / distance between
        # common SNPs), with the "distance between common SNPs" taken to be 1000 bp.
        d *= (self._read_length[sampidx] / 1000.)

      # Cap at 1e6 * read_depth.
      total_reads[sampidx] = int(round(min(d, self._max_reads(sampidx))))
    return total_reads

  def _format_overlapping_variants(self, variants, maj_cn, min_cn):
    assert len(set(maj_cn)) == len(set(min_cn)) == 1
    variants = [(ssm_id, str(min_cn[0]), str(maj_cn[0])) for ssm_id in variants]
    return variants

  def _format_cnvs(self, cnvs, variants):
    log('Estimated read depth: %s' % self._read_depth)

    for chrom, chrom_cnvs in cnvs.items():
      for cnv in chrom_cnvs:
        overlapping_variants = self._find_overlapping_variants(chrom, cnv, variants)
        total_reads = self._calc_total_reads(
          cnv['cell_prev'],
          cnv['start'],
          cnv['end'],
          cnv['major_cn'],
          cnv['minor_cn']
        )
        ref_reads = self._calc_ref_reads(cnv['cell_prev'], total_reads)
        yield {
          'chrom': chrom,
          'start': cnv['start'],
          'end': cnv['end'],
          'major_cn': cnv['major_cn'],
          'minor_cn': cnv['minor_cn'],
          'cellular_prevalence': cnv['cell_prev'],
          'ref_reads': ref_reads,
          'total_reads': total_reads,
          'overlapping_variants': self._format_overlapping_variants(overlapping_variants, cnv['major_cn'], cnv['minor_cn'])
        }

  def _merge_variants(self, cnv1, cnv2):
    cnv1_variant_names = set([v[0] for v in cnv1['overlapping_variants']])
    for variant in cnv2['overlapping_variants']:
      variant_name = variant[0]
      if variant_name not in cnv1_variant_names:
        cnv1['overlapping_variants'].append(variant)
      else:
        # If variant already in cnv1's list, ignore it. This should only occur
        # if two subclonal CNVs have close to 0.5 frequency each. In this case,
        # we lose information about major/minor status of the cnv2 relative to
        # its SSMs.
        log('%s already in %s' % (variant, cnv1['cnv_id']))

  # CNVs with similar a/d values should not be free to move around the
  # phylogeny independently, and so we merge them into a single entity. We may
  # do the same with SNVs bearing similar frequencies later on.
  def format_and_merge_cnvs(self, cnvs, variants, cellularity):
    formatted = list(self._format_cnvs(cnvs, variants))
    formatted.sort(key = lambda f: f['cellular_prevalence'][0], reverse = True)
    if len(formatted) == 0:
      return []

    for cnv in formatted:
      physical_cnvs = OrderedDict()
      for K in ('chrom', 'start', 'end', 'major_cn', 'minor_cn'):
        physical_cnvs[K] = cnv[K]

      assert len(set(physical_cnvs['major_cn'])) == len(set(physical_cnvs['major_cn'])) == 1
      physical_cnvs['major_cn'] = physical_cnvs['major_cn'][0]
      physical_cnvs['minor_cn'] = physical_cnvs['minor_cn'][0]

      physical_cnvs['cell_prev'] = '|'.join([str(C) for C in cnv['cellular_prevalence']])
      cnv['physical_cnvs'] = ','.join(['%s=%s' % (K, physical_cnvs[K]) for K in physical_cnvs.keys()])

    merged, formatted = formatted[:1], formatted[1:]
    merged[0]['cnv_id'] = 'c0'
    counter = 1

    for current in formatted:
      last = merged[-1]
      assert np.all(current['cellular_prevalence'] <= cellularity) and np.all(last['cellular_prevalence'] <= cellularity)

      # Only merge CNVs if they're clonal. If they're subclonal, leave them
      # free to move around the tree.
      if np.array_equal(current['cellular_prevalence'], last['cellular_prevalence']) \
      and np.array_equal(last['cellular_prevalence'], cellularity):
        # Merge the CNVs.
        log('Merging %s_%s and %s_%s' % (current['chrom'], current['start'], last['chrom'], last['start']))
        last['total_reads'] = current['total_reads'] + last['total_reads']
        last['ref_reads'] = self._calc_ref_reads(last['cellular_prevalence'], last['total_reads'])
        last['physical_cnvs'] += ';' + current['physical_cnvs']
        self._merge_variants(last, current)
      else:
        # Do not merge the CNVs.
        current['cnv_id'] = 'c%s' % counter
        merged.append(current)
        counter += 1

    for cnv in merged:
      cnv['ref_reads'] = np.round(cnv['ref_reads'] * self._cnv_confidence).astype(np.int)
      cnv['total_reads'] = np.round(cnv['total_reads'] * self._cnv_confidence).astype(np.int)

    return merged

class VariantFormatter(object):
  def __init__(self):
    self._counter = 0

  def _split_types(self, genotype):
    types = [int(e) for e in genotype.split('/')]
    if len(types) != 2:
      raise Exception('Not diploid: %s' % types)
    return types

  def _calc_ref_freq(self, ref_genotype, error_rate):
    types = self._split_types(ref_genotype)
    num_ref = len([t for t in types if t == 0])
    freq = (num_ref / 2) - error_rate
    if freq < 0:
      freq = 0.0
    if freq > 1:
      raise Exception('Nonsensical frequency: %s' % freq)
    return freq

  def format_variants(self, variants, ref_read_counts, total_read_counts, error_rate, sex):
    for variant_idx, variant in enumerate(variants):
      ssm_id = 's%s' % self._counter
      if hasattr(variant, 'ID') and variant.ID is not None:
        # This field will be defined by PyVCF, but not by our VariantId named
        # tuple that we have switched to, so this code will never actually run.
        # TODO: fix that.
        variant_name = variant.ID
      else:
        variant_name = '%s_%s' % (variant.CHROM, variant.POS)

      # TODO: switch back to using calc_ref_freq() when we no longer want mu_r
      # and mu_v fixed.
      # This is mu_r in PhyloWGS.
      expected_ref_freq = 1 - error_rate
      if variant.CHROM in ('Y', 'M') or (variant.CHROM == 'X' and sex == 'male'):
        # Haploid, so should only see non-variants when sequencing error
        # occurred. Note that chrY and chrM are always haploid; chrX is haploid
        # only in men, so script must know sex of patient to choose correct
        # value. Currently, I just assume that all data comes from men.
        #
        # This is mu_v in PhyloWGS.
        expected_var_freq = error_rate
      else:
        # Diploid, so should see variants in (0.5 - error_rate) proportion of
        # reads.
        #
        # This is mu_v in PhyloWGS.
        expected_var_freq = 0.5 - error_rate

      yield {
        'ssm_id': ssm_id,
        'chrom': variant.CHROM,
        'pos': variant.POS,
        'variant_name': variant_name,
        'ref_reads': list(ref_read_counts[variant_idx,:]),
        'total_reads': list(total_read_counts[variant_idx,:]),
        'expected_ref_freq': expected_ref_freq,
        'expected_var_freq': expected_var_freq,
      }
      self._counter += 1

def restricted_float(x):
  x = float(x)
  if x < 0.0 or x > 1.0:
    raise argparse.ArgumentTypeError('%r not in range [0.0, 1.0]' % x)
  return x

def chrom_key(chrom):
  if chrom.isdigit():
    return int(chrom)
  elif chrom == 'X':
    return 100
  elif chrom == 'Y':
    return 101
  else:
    raise Exception('Unknown chrom: %s' % chrom)

def variant_key(var):
  chrom = chrom_key(var.CHROM)
  return (chrom, var.POS)

class Segmenter(object):
  def _organize_cnvs(self, cnv_set):
    organized = defaultdict(list)

    for sampidx, cnvs in enumerate(cnv_set):
      for chrom, chrom_cnvs in cnvs.items():
        for cnv in chrom_cnvs:
          organized[chrom].append({
            'sample': sampidx,
            'start': cnv['start'],
            'end': cnv['end'],
            'major_cn': cnv['major_cn'],
            'minor_cn': cnv['minor_cn'],
            'cell_prev': cnv['cellular_prevalence']
          })

    for chrom, cnvs in organized.items():
      # Intervals may not be sorted in input file.
      cnvs.sort(key = lambda c: c['start'])

    return organized

  def _create_intervals(self, cnv_set):
    # intervals[chrom][(major, minor)]
    intervals = defaultdict(list)
    min_size_for_inclusion = 1

    for chrom, cnvs in cnv_set.items():
      for cnv in cnvs:
        # We sorted above to place start coordinates after end coordinates. But
        # if a CNV was listed with *the same* start and end position (meaning a
        # zero-length record, assuming intervals that are left-closed but
        # right-open), we will encounter the end for that record before its
        # start. As such, the "open_samples.remove()" call below will fail, as
        # the given intervals will not have been opened when we encounter its
        # end.
        #
        # Note the above assumes a half-open interpretation of intervals. I
        # don't think I implemented this -- if I recall, the code dealing with
        # CNVs (such as determining SSM overlap) assumes fully-closed intervals
        # (i.e., it doesn't check if cnv.start <= ssm.locus <= (cnv.end + 1)).
        # Normally this doesn't matter, given the low resolution of CNV calls
        # -- we should never encounter such small intervals. But a pathological
        # case in which CNV inputs had the same start & end coordinates for
        # some intervals revealed that the code crashes on this input. We
        # should provide a more informative error in such cases, which the
        # following assertion does.
        assert cnv['start'] < cnv['end'], ('In CNV %s, start position occurs at or after the end position' % cnv)

      start_pos = [(c['start'], 'start', (c['sample'], c['cell_prev'], c['major_cn'], c['minor_cn'])) for c in cnvs]
      end_pos   = [(c['end'],   'end',   (c['sample'], c['cell_prev'], c['major_cn'], c['minor_cn'])) for c in cnvs]

      # True > False, so this sorting will place start positions after end
      # positions if both have same coordinate.
      positions = sorted(start_pos + end_pos, key = lambda e: (e[0], e[1] == 'start'))
      assert len(positions) >= 2, 'Fewer than two positions in %s' % positions

      # prev_pos is updated each time we move to a new coordinate on the
      # chromosome. Multiple start or end points may be associated with any
      # given coordinate.
      prev_pos = None
      open_samples = []
      idx = 0

      while idx < len(positions):
        points_at_locus = [positions[idx]]
        locus = points_at_locus[0][0]

        # Gather all interval breakpoints at this locus.
        while True:
          assert positions[idx][0] >= locus
          idx += 1
          if idx == len(positions) or positions[idx][0] > locus:
            break
          points_at_locus.append(positions[idx])

        if prev_pos is None:
          assert len(open_samples) == 0

        if len(open_samples) > 0:
          # If some samples are already open from previous loci (such that
          # last_pos will not be None), add this interval.
          assert locus > prev_pos
          interval = (prev_pos, locus)
          if interval[1] - interval[0] > min_size_for_inclusion:
            intervals[chrom].append((interval[0], interval[1], sorted(open_samples)))
        else:
          # All points should be start points.
          assert set([i[1] for i in points_at_locus]) == set(['start'])

        prev_pos = locus

        # Update open_samples in accordance with whether each breakpoint at
        # this locus starts or ends an interval.
        for pos, pt_type, (sampidx, cell_prev, major_cn, minor_cn) in points_at_locus:
          if pt_type == 'start':
            log('Adding ', (pos, pt_type, sampidx, cell_prev, major_cn, minor_cn))
            open_samples.append((sampidx, cell_prev, major_cn, minor_cn))
          elif pt_type == 'end':
            log('Removing ', (pos, pt_type, sampidx, cell_prev, major_cn, minor_cn))
            open_samples.remove((sampidx, cell_prev, major_cn, minor_cn))
          else:
            raise Exception('Unknown point type: %s' % pt_type)

      assert len(open_samples) == 0

    return intervals

  def _merge_adjacent(self, cncalls, allowed_gap = 0):
    cncalls.sort(key = lambda c: (Util.chrom_key(c['chrom']), c['start']))
    merged = []
    idx = 0
    while idx < len(cncalls):
      adjacent = [cncalls[idx]]
      idx += 1

      while idx < len(cncalls) and \
      cncalls[idx]['chrom'] == adjacent[-1]['chrom'] and \
      cncalls[idx]['major_cn'] == adjacent[-1]['major_cn'] and \
      cncalls[idx]['minor_cn'] == adjacent[-1]['minor_cn'] and \
      0 <= cncalls[idx]['start'] - adjacent[-1]['end'] <= allowed_gap:
        adjacent.append(cncalls[idx])
        idx += 1

      if len(adjacent) > 1:
        log('Merging ', adjacent)
        copy = dict(adjacent[0])
        copy['end'] = adjacent[-1]['end']
        merged.append(copy)
      else:
        merged.append(adjacent[0])

    return merged

  def segment(self, cn_calls):
    # Merge adjacent CNVs here rather than when data loaded, as what can be
    # merged will be determined by what tetraploidy correction, if any, is
    # applied to the data.
    #for sampidx, cnvs in enumerate(cn_calls):
      #cn_calls[sampidx] = self._merge_adjacent(cnvs)
    organized = self._organize_cnvs(cn_calls)
    return self._create_intervals(organized)

class MultisampleCnvCombiner(object):
  def __init__(self, cn_regions, cellularity):
    self.sampidxs = set(range(len(cn_regions)))
    segments = Segmenter().segment(cn_regions)
    self._cnvs = self._reformat_segments_as_cnvs(segments)
    self._cellularity = cellularity

  def _reformat_segments_as_cnvs(self, segments):
    reformatted = defaultdict(list)
    _retrieve_val = lambda idx: np.array(zip(*open_samples)[idx])

    for chrom, chrom_cnvs in segments.items():
        for start, end, open_samples in chrom_cnvs:
          sampidx = _retrieve_val(0)
          cell_prev = _retrieve_val(1)
          major_cn = _retrieve_val(2)
          minor_cn = _retrieve_val(3)
          cnv = {
            'start': start,
            'end': end,
            'cell_prev': cell_prev,
            'major_cn': major_cn,
            'minor_cn': minor_cn,
            'sampidx': sampidx,
          }
          reformatted[chrom].append(cnv)

    return reformatted

  def _ensure_no_overlap(self, cnvs):
    for chrom, chrom_cnvs in cnvs.items():
      for idx in range(len(chrom_cnvs) - 1):
        current, next = chrom_cnvs[idx], chrom_cnvs[idx + 1]
        assert current['start'] < current['end'] <= next['start'] < next['end']

  def _is_region_normal_cn(self, major, minor):
    return major == minor == 1

  def _is_multisample_region_normal_cn(self, major, minor):
    return set(major) == set(minor) == set([1])


  def _get_abnormal_state_for_all_samples(self, cnv):
    # All samples must have at least one record for this region, or don't
    # include it.
    if set(cnv['sampidx']) != self.sampidxs:
      return None

    abnormal_state = None
    filtered = []

    for sampidx, cell_prev, major, minor in zip(cnv['sampidx'], cnv['cell_prev'], cnv['major_cn'], cnv['minor_cn']):
      # Region may be (clonal or subclonal) normal in a sample, so ignore such records.
      if self._is_region_normal_cn(major, minor):
        continue

      # Either we haven't observed an abnormal CN state in this region before,
      # or the observed abnormal state matches what we've already seen.
      if abnormal_state is None or abnormal_state == (major, minor):
        abnormal_state = (major, minor)
        filtered.append({'sampidx': sampidx, 'cell_prev': cell_prev, 'major_cn': major, 'minor_cn': minor})
        continue
      # The abnormal state (i.e., major & minor alleles) is *different* from
      # what we've seen before. The PWGS model doesn't currently account for
      # such cases, so ignore the region.
      else:
        return None

    # None of the observed records were abnormal -- i.e., all samples report
    # the region is normal. Reject the region.
    if abnormal_state is None:
      return None

    retained_sampidxs = [F['sampidx'] for F in filtered]
    # Sanity check: when we originally parsed the CNVs, the samples should have
    # been added in order, and that ought not to have changed.
    assert retained_sampidxs == sorted(retained_sampidxs)
    # Sanity check: we should have no duplicate samples. While a given sample
    # may report any number of records for a region, above we discarded normal
    # regions, and ensured that only one abnormal state exists in all samples.
    # Thus, we should have no more than one record per sample for this region.
    assert len(retained_sampidxs) == len(set(retained_sampidxs))

    # Add a record for all samples that reported this region as clonal normal.
    cell_prev_when_absent = 0
    for missing_sampidx in self.sampidxs - set(retained_sampidxs):
      filtered.append({
        'sampidx': missing_sampidx,
        'cell_prev': cell_prev_when_absent,
        'major_cn': abnormal_state[0],
        'minor_cn': abnormal_state[1]
      })
    # Sort by sampidx.
    filtered.sort(key = lambda F: F['sampidx'])
    # Ensure all samples have one record.
    assert len(filtered) == len(self.sampidxs)

    return filtered

  def load_single_abnormal_state_cnvs(self):
    '''
    Return all regions that possess at most one abnormal state across samples.
    E.g., given three samples, S_1 and S_3 report the region as (2, 1) (with
    potentially different cellular prevalences), while S_2 lists it as clonal
    (1, 1). In such an instance, the record for S_2 will *not* indicate the
    region is normal. Instead, the S_2 record will show a state of (2, 1) with
    a cellular prevalence of zero. This is done so that we can calculate
    sensible `a` and `d` values for cnv_data.txt.
    '''
    # In Battenberg, either one region is normal and the other abnormal,
    # or both are abnormal.
    # In TITAN, only one abnormal region will be listed, without a
    # corresponding normal region.
    abnormal_cnvs = defaultdict(list)

    for chrom, chrom_cnvs in self._cnvs.items():
      if not is_good_chrom(chrom):
        continue
      for cnv in chrom_cnvs:
        states_for_all_samples = self._get_abnormal_state_for_all_samples(cnv)
        if states_for_all_samples is None:
          continue

        combined_states = { K: np.array([S[K] for S in states_for_all_samples]) for K in states_for_all_samples[0].keys() }
        cnv.update(combined_states)
        abnormal_cnvs[chrom].append(cnv)
      abnormal_cnvs[chrom].sort(key = lambda C: C['start'])

    self._ensure_no_overlap(abnormal_cnvs)
    return abnormal_cnvs

  def load_normal_cnvs(self):
    '''
    Return all regions that are clonal normal across all samples.
    '''
    normal_cnvs = defaultdict(list)

    for chrom, chrom_cnvs in self._cnvs.items():
      if not is_good_chrom(chrom):
        continue
      for cnv in chrom_cnvs:
        if not self._is_multisample_region_normal_cn(cnv['major_cn'], cnv['minor_cn']):
          continue
        if not set(cnv['sampidx']) == self.sampidxs:
          continue
        if not np.array_equal(cnv['cell_prev'], self._cellularity):
          # The region must be clonal normal to be retained. This check
          # shouldn't be necessary, as we've already ensured all calls have
          # major = minor = 1, but we perform it just to be thorough.
          continue
        normal_cnvs[chrom].append(cnv)
      normal_cnvs[chrom].sort(key = lambda C: C['start'])

    self._ensure_no_overlap(normal_cnvs)
    return normal_cnvs

  def load_cnvs(self):
    '''
    Return all regions that are clonal normal across all samples, or possess at
    most one (clonal or subclonal) abnormal state in each sample.
    '''
    combined = defaultdict(list)

    normal_cnvs = self.load_normal_cnvs()
    abnormal_cnvs = self.load_single_abnormal_state_cnvs()
    for chrom in set(normal_cnvs.keys()) | set(abnormal_cnvs.keys()):
      combined[chrom] = normal_cnvs[chrom] + abnormal_cnvs[chrom]
      combined[chrom].sort(key = lambda C: C['start'])
    self._ensure_no_overlap(combined)

    return combined

class VariantAndCnvGroup(object):
  def __init__(self):
    self._multisamp_cnv = None
    self._cellularity = None

  def add_variants(self, variants, ref_read_counts, total_read_counts):
    self._variants = variants
    self._variant_idxs = list(range(len(variants)))
    self._ref_read_counts = ref_read_counts
    self._total_read_counts = total_read_counts
    # Estimate read depth before any filtering of variants is performed, in
    # case no SSMs remain afterward.
    self._estimated_read_depth = self._estimate_read_depth()

  def _find_cellularity(self, cnvs):
    max_cellular_prevs = np.zeros(len(cnvs))

    for sampidx, sample_cnvs in enumerate(cnvs):
      for chrom_regions in sample_cnvs.values():
        for cnr in chrom_regions:
          if cnr['cellular_prevalence'] > max_cellular_prevs[sampidx]:
            max_cellular_prevs[sampidx] = cnr['cellular_prevalence']

    return max_cellular_prevs

  def add_cnvs(self, cn_regions):
    self._cellularity = self._find_cellularity(cn_regions)
    self._multisamp_cnv = MultisampleCnvCombiner(cn_regions, self._cellularity)
    self._sampidxs = self._multisamp_cnv.sampidxs

  def has_cnvs(self):
    return self._multisamp_cnv is not None

  def _filter_variants_outside_regions(self, regions, before_label, after_label):
    def _is_pos_in_regions(chrom, pos):
      for cnv in regions[chrom]:
        if cnv['start'] <= pos <= cnv['end']:
          return True
      return False

    filtered = []

    for vidx in self._variant_idxs:
      variant = self._variants[vidx]
      if _is_pos_in_regions(variant.CHROM, variant.POS):
        filtered.append(vidx)

    self._print_variant_differences(
      [self._variants[idx] for idx in self._variant_idxs],
      [self._variants[idx] for idx in filtered],
      before_label,
      after_label
    )
    self._variant_idxs = filtered

  def _print_variant_differences(self, before, after, before_label, after_label):
    before = set(before)
    after = set(after)
    log('%s=%s %s=%s delta=%s' % (before_label, len(before), after_label, len(after), len(before) - len(after)))

    assert after.issubset(before)
    removed = list(before - after)
    removed.sort(key = variant_key)

    def _print_region(var):
      var_name = '%s_%s' % (var.CHROM, var.POS)
      region_type = None
      containing_cnv = None

      for cnv in self._multisamp_cnv.load_normal_cnvs()[var.CHROM]:
        if cnv['start'] <= var.POS <= cnv['end']:
          region_type = 'normal'
          containing_cnv = cnv
          break
      for cnv in self._multisamp_cnv.load_single_abnormal_state_cnvs()[var.CHROM]:
        if cnv['start'] <= var.POS <= cnv['end']:
          assert region_type is None and containing_cnv is None
          region_type = 'abnormal'
          containing_cnv = cnv
          break

      if containing_cnv is not None:
        log('%s\t[in %s-CN region chr%s(%s, %s)]' % (var_name, region_type, var.CHROM, containing_cnv['start'], containing_cnv['end']))
      else:
        log('%s\t[outside all regions]' % var_name)

    for var in removed:
      _print_region(var)

  def retain_only_variants_in_normal_cn_regions(self):
    if not self.has_cnvs():
      raise Exception('CN regions not yet provided')

    normal_cn = self._multisamp_cnv.load_normal_cnvs()
    filtered = self._filter_variants_outside_regions(normal_cn, 'all_variants', 'only_normal_cn')

  def exclude_variants_in_multiple_abnormal_or_unlisted_regions(self):
    # Battenberg:
    #   Five possible placements for variant in Battenberg according to CN records:
    #   1 record:
    #     That record has normal CN: include
    #     That record has abnormal CN: include
    #   2 records:
    #     One record is normal CN, one record is abnormal CN: include
    #     Both records are abnormal CN: exclude (as we don't know what order the CN events occurred in)
    # TITAN:
    #   In output seen to date, TITAN will only list one record per region. If
    #   the CN state is abnormal and clonal_frac < 1, this implies the
    #   remainder of the region will be normal CN. Multiple abnormal records
    #   for the same region are likely possible, but I haven't yet seen any.
    #   Regardless, when they occur, they should be properly handled by the
    #   code.
    if not self.has_cnvs():
      raise Exception('CN regions not yet provided')

    # If variant isn't listed in *any* region: exclude (as we suspect CNV
    # caller didn't know what to do with the region).
    self._filter_variants_outside_regions(self._multisamp_cnv.load_cnvs(), 'all_variants', 'outside_subclonal_cn')

  def format_variants(self, sample_size, error_rate, priority_ssms, sex):
    if sample_size is None:
      sample_size = len(self._variant_idxs)
    random.shuffle(self._variant_idxs)

    subsampled, remaining = [], []

    for variant_idx in self._variant_idxs:
      variant = self._variants[variant_idx]
      if len(subsampled) < sample_size and (variant.CHROM, variant.POS) in priority_ssms:
        subsampled.append(variant_idx)
      else:
        remaining.append(variant_idx)

    assert len(subsampled) <= sample_size
    needed = sample_size - len(subsampled)
    subsampled = subsampled + remaining[:needed]
    nonsubsampled = remaining[needed:]

    subsampled.sort(key = lambda idx: variant_key(self._variants[idx]))
    subsampled_variants = get_elements_at_indices(self._variants, subsampled)
    subsampled_ref_counts = self._ref_read_counts[subsampled,:]
    subsampled_total_counts = self._total_read_counts[subsampled,:]

    nonsubsampled.sort(key = lambda idx: variant_key(self._variants[idx]))
    nonsubsampled_variants = get_elements_at_indices(self._variants, nonsubsampled)
    nonsubsampled_ref_counts = self._ref_read_counts[nonsubsampled,:]
    nonsubsampled_total_counts = self._total_read_counts[nonsubsampled,:]

    formatter = VariantFormatter()
    subsampled_formatted = list(formatter.format_variants(subsampled_variants, subsampled_ref_counts, subsampled_total_counts, error_rate, sex))
    nonsubsampled_formatted = list(formatter.format_variants(nonsubsampled_variants, nonsubsampled_ref_counts, nonsubsampled_total_counts, error_rate, sex))

    return (subsampled_formatted, nonsubsampled_formatted)

  def write_variants(self, variants, outfn):
    with open(outfn, 'w') as outf:
      print('\t'.join(('id', 'gene', 'a', 'd', 'mu_r', 'mu_v')), file=outf)
      for variant in variants:
        variant['ref_reads'] = ','.join([str(v) for v in variant['ref_reads']])
        variant['total_reads'] = ','.join([str(v) for v in variant['total_reads']])
        vals = (
          'ssm_id',
          'variant_name',
          'ref_reads',
          'total_reads',
          'expected_ref_freq',
          'expected_var_freq',
        )
        vals = [variant[k] for k in vals]
        print('\t'.join([str(v) for v in vals]), file=outf)

  def _estimate_read_depth(self):
    read_sum = 0
    if len(self._variants) == 0:
      default_read_depth = 50
      log('No variants available, so fixing read depth at %s.' % default_read_depth)
      return default_read_depth
    else:
      return np.nanmean(self._total_read_counts, axis=0)

  def write_cnvs(self, variants, outfn, cnv_confidence, read_length):
    with open(outfn, 'w') as outf:
      print('\t'.join(('cnv', 'a', 'd', 'ssms', 'physical_cnvs')), file=outf)
      formatter = CnvFormatter(cnv_confidence, self._estimated_read_depth, read_length, self._sampidxs)
      # Last place I'm working is at this position
      for cnv in formatter.format_and_merge_cnvs(self._multisamp_cnv.load_single_abnormal_state_cnvs(), variants, self._cellularity):
        overlapping = [','.join(o) for o in cnv['overlapping_variants']]
        vals = (
          cnv['cnv_id'],
          ','.join([str(V) for V in cnv['ref_reads']]),
          ','.join([str(V) for V in cnv['total_reads']]),
          ';'.join(overlapping),
          cnv['physical_cnvs']
        )
        print('\t'.join(vals), file=outf)

def log(*msgs):
  if log.verbose:
    print(*msgs, file=sys.stderr)
log.verbose = False

class CnvParser(object):
  def __init__(self, cn_filename):
    self._cn_filename = cn_filename

  def parse(self):
    cn_regions = defaultdict(list)

    with open(self._cn_filename) as cnf:
      reader = csv.DictReader(cnf, delimiter='\t')
      for record in reader:
        chrom = record['chromosome'].upper()
        del record['chromosome']
        for key in ('start', 'end', 'major_cn', 'minor_cn'):
          # Some records from Battenberg have major and minor listed as, e.g.,
          # "1.0", so cast to float before int.
          assert float(record[key]) == int(float(record[key]))
          record[key] = int(float(record[key]))
        record['cellular_prevalence'] = float(record['cellular_prevalence'])
        cn_regions[chrom].append(record)

    # Ensure CN regions are properly sorted, which we later rely on when
    # filtering out regions with multiple abnormal CN states.
    for chrom, regions in cn_regions.items():
      cn_regions[chrom] = sorted(regions, key = lambda r: r['start'])

    return cn_regions

def get_elements_at_indices(L, indices):
  elem = []
  for idx in indices:
    elem.append(L[idx])
  return elem

def parse_priority_ssms(priority_ssm_filename):
  if priority_ssm_filename is None:
    return set()
  with open(priority_ssm_filename) as priof:
    lines = priof.readlines()

  priority_ssms = set()
  for line in lines:
    chrom, pos = line.strip().split('_', 1)
    priority_ssms.add((chrom.upper(), int(pos)))
  return set(priority_ssms)

def impute_missing_total_reads(total_reads, missing_variant_confidence):
  # Change NaNs to masked values via SciPy.
  masked_total_reads = ma.fix_invalid(total_reads)

  # Going forward, suppose you have v variants and s samples in a v*s matrix of
  # read counts. Missing values are masked.

  # Calculate geometric mean of variant read depth in each sample. Result: s*1
  sample_means = gmean(masked_total_reads, axis=0)
  assert np.sum(sample_means <= 0) == np.sum(np.isnan(sample_means)) == 0
  # Divide every variant's read count by its mean sample read depth to get read
  # depth enrichment relative to other variants in sample. Result: v*s
  normalized_to_sample = np.dot(masked_total_reads, np.diag(1./sample_means))
  # For each variant, calculate geometric mean of its read depth enrichment
  # across samples. Result: v*1
  variant_mean_reads = gmean(normalized_to_sample, axis=1)
  assert np.sum(variant_mean_reads <= 0) == np.sum(np.isnan(variant_mean_reads)) == 0

  # Convert 1D arrays to vectors to permit matrix multiplication.
  imputed_counts = np.dot(variant_mean_reads.reshape((-1, 1)), sample_means.reshape((1, -1)))
  nan_coords = np.where(np.isnan(total_reads))
  total_reads[nan_coords] = imputed_counts[nan_coords]
  assert np.sum(total_reads <= 0) == np.sum(np.isnan(total_reads)) == 0

  total_reads[nan_coords] *= missing_variant_confidence
  return np.floor(total_reads).astype(np.int)

def impute_missing_ref_reads(ref_reads, total_reads):
  ref_reads = np.copy(ref_reads)

  assert np.sum(np.isnan(total_reads)) == 0
  nan_coords = np.where(np.isnan(ref_reads))
  ref_reads[nan_coords] = total_reads[nan_coords]
  assert np.sum(np.isnan(ref_reads)) == 0

  return ref_reads.astype(np.int)

def is_good_chrom(chrom):
  # Ignore the following:
  #   * Variants unmapped ('chrUn') or mapped to fragmented chromosome ('_random')
  #   * Weird chromosomes from Mutect (e.g., "chr17_ctg5_hap1").
  #   * Mitochondrial ("mt" or "m"), which are weird
  #   * Sex chromosomes difficult to deal with, as expected frequency depends on
  #     whether patient is male or female, so ignore them for now. TODO: fix this.
  if chrom in [str(i) for i in range(1, 23)] + ['X', 'Y']:
    return True
  else:
    return False

def parse_variants(args, vcf_types, num_samples):
  parsed_variants = []
  all_variant_ids = []

  VariantId = namedtuple('VariantId', ['CHROM', 'POS'])

  for vcf_file in args.vcf_files:
    vcf_type, vcf_fn = vcf_file.split('=', 1)

    if vcf_type == 'sanger':
      variant_parser = SangerParser(vcf_fn, args.tumor_sample)
    elif vcf_type == 'mutect_pcawg':
      variant_parser = MutectPcawgParser(vcf_fn, args.tumor_sample)
    elif vcf_type == 'mutect_smchet':
      variant_parser = MutectSmchetParser(vcf_fn, args.tumor_sample)
    elif vcf_type == 'mutect_tcga':
      variant_parser = MutectTcgaParser(vcf_fn, args.tumor_sample)
    elif vcf_type == 'muse':
      variant_parser = MuseParser(vcf_fn, args.muse_tier, args.tumor_sample)
    elif vcf_type == 'dkfz':
      variant_parser = DKFZParser(vcf_fn, args.tumor_sample)
    elif vcf_type == 'strelka':
      variant_parser = StrelkaParser(vcf_fn, args.tumor_sample)
    elif vcf_type == 'vardict':
      variant_parser = VarDictParser(vcf_fn, args.tumor_sample)
    elif vcf_type == 'pcawg_consensus':
      variant_parser = PcawgConsensusParser(vcf_fn, args.tumor_sample)
    else:
      raise Exception('Unknowon variant type: %s' % vcf_type)

    parsed_variants.append(variant_parser.list_variants())
    variant_ids = [VariantId(str(v[0].CHROM), int(v[0].POS)) for v in parsed_variants[-1]]
    all_variant_ids += variant_ids

  all_variant_ids = list(set(all_variant_ids)) # Eliminate duplicates.
  all_variant_ids.sort(key = variant_key)
  num_variants = len(all_variant_ids)
  variant_positions = dict(zip(all_variant_ids, range(num_variants)))

  total_read_counts = np.zeros((num_variants, num_samples))
  total_read_counts.fill(np.nan)
  ref_read_counts = np.copy(total_read_counts)

  for sample_idx, parsed in enumerate(parsed_variants):
    for variant, ref_reads, total_reads in parsed:
      variant_id = VariantId(str(variant.CHROM), int(variant.POS))
      variant_idx = variant_positions[variant_id]
      ref_read_counts[variant_idx, sample_idx] = ref_reads
      total_read_counts[variant_idx, sample_idx] = total_reads

  total_read_counts = impute_missing_total_reads(total_read_counts, args.missing_variant_confidence)
  ref_read_counts = impute_missing_ref_reads(ref_read_counts, total_read_counts)
  return (all_variant_ids, ref_read_counts, total_read_counts)

def infer_sex(variant_ids):
  num_y_variants = len([V for V in variant_ids if V.CHROM == 'Y'])
  if num_y_variants > 0:
    return 'male'
  else:
    return 'female'

def main():
  vcf_types = set(('sanger', 'mutect_pcawg', 'mutect_smchet', 'mutect_tcga', 'muse','dkfz', 'strelka', 'vardict', 'pcawg_consensus'))

  parser = argparse.ArgumentParser(
    description='Create ssm_dat.txt and cnv_data.txt input files for PhyloWGS from VCF and CNV data.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('-e', '--error-rate', dest='error_rate', type=restricted_float, default=0.001,
    help='Expected error rate of sequencing platform')
  parser.add_argument('--missing-variant-confidence', dest='missing_variant_confidence', type=restricted_float, default=1.,
    help='Confidence in range [0, 1] that SSMs missing from a sample are indeed not present in that sample')
  parser.add_argument('-s', '--sample-size', dest='sample_size', type=int,
    help='Subsample SSMs to reduce PhyloWGS runtime')
  parser.add_argument('-P', '--priority-ssms', dest='priority_ssm_filename',
    help='File containing newline-separated list of SSMs in "<chr>_<locus>" format to prioritize for inclusion')
  parser.add_argument('--cnvs', dest='cnv_files', action='append',
    help='Path to per-sample CNV files created with parse_cnvs.py. Specified once per CNV file.')
  parser.add_argument('--regions', dest='regions', choices=('normal_cn', 'normal_and_abnormal_cn', 'all'), default='normal_and_abnormal_cn',
    help='Which regions to use variants from. Refer to the parser README for more details.')
  parser.add_argument('--output-cnvs', dest='output_cnvs', default='cnv_data.txt',
    help='Output destination for CNVs')
  parser.add_argument('--output-variants', dest='output_variants', default='ssm_data.txt',
    help='Output destination for variants')
  parser.add_argument('--tumor-sample', dest='tumor_sample',
    help='Name of the tumor sample in the input VCF file. Defaults to last sample if not specified.')
  parser.add_argument('--cnv-confidence', dest='cnv_confidence', type=restricted_float, default=0.5,
    help='Confidence in CNVs. Set to < 1 to scale "d" values used in CNV output file')
  parser.add_argument('--read-length', dest='read_length', type=int, default=100,
    help='Approximate length of reads. Used to calculate confidence in CNV frequencies')
  parser.add_argument('--muse-tier', dest='muse_tier', type=int, default=0,
    help='Maximum MuSE tier to include')
  parser.add_argument('--nonsubsampled-variants', dest='output_nonsubsampled_variants',
    help='If subsampling, write nonsubsampled variants to separate file, in addition to subsampled variants')
  parser.add_argument('--nonsubsampled-variants-cnvs', dest='output_nonsubsampled_variants_cnvs',
    help='If subsampling, write CNVs for nonsubsampled variants to separate file')
  parser.add_argument('--sex', dest='sex', default='auto', choices=('auto', 'male', 'female'),
    help='Sex of patient. Used to adjust expected variant frequencies on sex chromosomes. ' +
    'If auto, patient is set to male if any variants are provided on the Y chromosome, and female otherwise.')
  parser.add_argument('--verbose', dest='verbose', action='store_true')
  parser.add_argument('vcf_files', nargs='+', help='One or more space-separated occurrences of <vcf_type>=<path>. E.g., sanger=variants1.vcf muse=variants2.vcf. Valid vcf_type values: %s' % ', '.join(vcf_types))
  args = parser.parse_args()

  log.verbose = args.verbose
  num_samples = len(args.vcf_files)

  variant_ids, ref_read_counts, total_read_counts = parse_variants(args, vcf_types, num_samples)

  # Fix random seed to ensure same set of SSMs chosen when subsampling on each
  # invocation.
  random.seed(1)

  grouper = VariantAndCnvGroup()
  grouper.add_variants(variant_ids, ref_read_counts, total_read_counts)

  if args.cnv_files:
    if len(args.cnv_files) != len(args.vcf_files):
      raise Exception('Number of CNV samples must match number of VCF samples. You provided %s CNV sample(s) and %s VCF sample(s)' % (len(args.cnv_files), len(args.vcf_files)))
    cn_regions = [CnvParser(cnvf).parse() for cnvf in args.cnv_files]
    grouper.add_cnvs(cn_regions)

  if not grouper.has_cnvs():
    assert args.regions == 'all', 'If you do not provide CNA data, you must specify --regions=all'

  if args.regions == 'normal_cn':
    grouper.retain_only_variants_in_normal_cn_regions()
  elif args.regions == 'normal_and_abnormal_cn':
    grouper.exclude_variants_in_multiple_abnormal_or_unlisted_regions()
  elif args.regions == 'all':
    pass
  else:
    raise Exception('Unknown --regions value: %s' % args.regions)

  priority_ssms = parse_priority_ssms(args.priority_ssm_filename)

  if args.sex == 'auto':
    sex = infer_sex(variant_ids)
  else:
    sex = args.sex

  subsampled_vars, nonsubsampled_vars = grouper.format_variants(args.sample_size, args.error_rate, priority_ssms, sex)
  if len(subsampled_vars) == 0:
    raise Exception('No variants to write')
  grouper.write_variants(subsampled_vars, args.output_variants)
  if args.output_nonsubsampled_variants:
    grouper.write_variants(nonsubsampled_vars, args.output_nonsubsampled_variants)

  if grouper.has_cnvs() and args.regions != 'normal_cn':
    # Write CNVs.
    read_length = num_samples * [args.read_length]
    grouper.write_cnvs(subsampled_vars, args.output_cnvs, args.cnv_confidence, read_length)
    if args.output_nonsubsampled_variants and args.output_nonsubsampled_variants_cnvs:
      grouper.write_cnvs(nonsubsampled_vars, args.output_nonsubsampled_variants_cnvs, args.cnv_confidence, read_length)
  else:
    # Write empty CNV file.
    with open(args.output_cnvs, 'w'):
      pass

if __name__ == '__main__':
  main()
