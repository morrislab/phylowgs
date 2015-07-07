#!/usr/bin/env python2
from __future__ import print_function

# Requires PyVCF. To install: pip3 install pyvcf
import vcf
import argparse
from collections import namedtuple, defaultdict, OrderedDict
import random

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
      ref_reads, total_reads = self._calc_read_counts(variant)
      variants_and_reads.append((variant, ref_reads, total_reads))
    return variants_and_reads

  def _calc_read_counts(self, variant):
    raise Exception('Not implemented -- use child class')

  def _parse_vcf(self, vcf_filename):
    vcfr = vcf.Reader(filename=vcf_filename)
    records = []
    for variant in vcfr:
      variant.CHROM = variant.CHROM.lower()
      # Some VCF dialects prepend "chr", some don't. Remove the prefix to
      # standardize.
      if variant.CHROM.startswith('chr'):
        variant.CHROM = variant.CHROM[3:]
      records.append(variant)
    return records

  def _is_good_chrom(self, chrom):
    if chrom.startswith('un') or chrom.endswith('_random'):
      # Variant unmapped ('chrUn') or mapped to fragmented chromosome
      # ('_random'), so ignore.
      return False

    if chrom.startswith('hs37d5') or chrom.startswith('gl'):
      return False

    if '_' in chrom:
      # Reject weird chromosomes from Mutect (e.g., "chr17_ctg5_hap1").
      return False

    # Mitochondrial are weird, so ignore them.
    if chrom in ('m', 'mt'):
      return False
    # Sex chromosomes difficult to deal with, as expected frequency depends on
    # whether patient is male or female, so ignore them for now.
    #
    # TODO: properly deal with sex chromsomes.
    if chrom in ('x', 'y'):
      return False

    return True

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
      if not self._is_good_chrom(variant.CHROM):
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
      return -1

class SangerParser(VariantParser):
  '''
  Works with Battenberg-formatted CNV calls.
  '''
  def __init__(self, vcf_filename):
    self._vcf_filename = vcf_filename

  def _find_ref_and_variant_nt(self, variant):
    # Try most probable genotype called by CaVEMan. If can't find variant nt in
    # it, try the second most probable genotype.
    genotypes = [variant.INFO['TG'][0], variant.INFO['SG'][0]]
    variant_set = set()

    while len(variant_set) == 0:
      if len(genotypes) == 0:
        raise Exception('No more genotypes to find variant_nt in for %s' % variant)
      genotype = genotypes.pop(0)
      normal_genotype, tumor_genotype = genotype.split('/')
      # TODO: We throw out hetero germline. Deal with this later.
      if normal_genotype[0] != normal_genotype[1]:
        print('Ignoring heterozygous normal genotype %s' % normal_genotype, file=sys.stderr)
      reference_nt = normal_genotype[0]
      variant_set = set(tumor_genotype) - set(reference_nt)

    variant_nt = variant_set.pop()
    return (reference_nt, variant_nt)

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

class MuseParser(VariantParser):
  def __init__(self, vcf_filename, tier=0, tumor_sample=None):
    self._vcf_filename = vcf_filename
    self._tier = tier
    self._tumor_sample = tumor_sample

  def _calc_read_counts(self, variant):
    tumor_i = self._get_tumor_index(variant, self._tumor_sample)
    total_reads = int(variant.samples[tumor_i]['DP4'])
    variant_reads = int(variant.samples[tumor_i]['AD'])

    ref_reads = total_reads - variant_reads

    return (ref_reads, total_reads)

  def _does_variant_pass_filters(self, variant):
    if variant.FILTER is None or len(variant.FILTER) == 0:
      return True
    if int(variant.FILTER[0][-1]) <= self._tier:
      # Variant failed one or more filters.
      return True
    return False
    
class StrelkaParser(VariantParser):
  def __init__(self, vcf_filename):
    self._vcf_filename = vcf_filename
    
  def _does_variant_pass_filters(self, variant):
    # Strelka outputs two files one for SNPs, the other for InDels
    # For now only deal with SNP file from Strelka
    if variant.is_snp:
      if variant.FILTER is None or len(variant.FILTER) == 0: 
        return True
    return False

  def _calc_read_counts(self, variant):
    alt = variant.ALT[0]
    total_reads = int(variant.samples[-1]['DP'])
	
    if alt is None:
      total_reads = 0
      variant_reads = 0
    else:
      variant_reads = getattr(variant.samples[-1].data, str(alt)+'U')[0]

    ref_reads = total_reads - variant_reads
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
    ref_reads = variant.samples[tumor_i]['AD'][0]
    variant_reads = variant.samples[tumor_i]['AD'][1]
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
    tumor_i = self._get_tumor_index(variant, self._tumor_sample)
    for_ref_reads = variant.samples[tumor_i]['DP4'][0]
    back_ref_reads = variant.samples[tumor_i]['DP4'][1]
    for_variant_reads = variant.samples[tumor_i]['DP4'][2]
    back_variant_reads = variant.samples[tumor_i]['DP4'][3]
    ref_reads = for_ref_reads + back_ref_reads
    var_reads = for_variant_reads + back_variant_reads
    total_reads = ref_reads + var_reads

    return (ref_reads, total_reads)

class BattenbergParser(object):
  def __init__(self, bb_filename):
    self._bb_filename = bb_filename

  def _compute_cn(self, cnv1, cnv2):
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
        chrom = fields[1].lower()
        start = int(fields[2])
        end = int(fields[3])
        pval = float(fields[5])

        cnv1 = OrderedDict()
        cnv1['start'] = start
        cnv1['end'] = end
        cnv1['nmaj'] = int(fields[8])
        cnv1['nmin'] = int(fields[9])
        cnv1['frac'] = float(fields[10])

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
          cnv2 = OrderedDict()
          cnv2['start'] = start
          cnv2['end'] = end
          cnv2['nmaj'] = int(fields[11])
          cnv2['nmin'] = int(fields[12])
          cnv2['frac'] = float(fields[13])
        else:
          cnv1['frac'] = 1.0

        cn_regions[chrom].append(cnv1)
        if cnv2 is not None:
          cn_regions[chrom].append(cnv2)
    return cn_regions

class CnvFormatter(object):
  def __init__(self, cnv_confidence, cellularity, read_depth, read_length):
    self._cnv_confidence = cnv_confidence
    self._cellularity = cellularity
    self._read_depth = read_depth
    self._read_length = read_length

  def _max_reads(self):
    return 1e6 * self._read_depth

  def _find_overlapping_variants(self, chrom, cnv, variants):
    overlapping = []

    start = cnv['start']
    end = cnv['end']
    for variant in variants:
      if chrom.lower() == variant['chrom'].lower():
        if start <= variant['pos'] <= end:
          overlapping.append(variant['ssm_id'])
    return overlapping

  def _calc_ref_reads(self, pop_frac, total_reads):
    tumor_cells_frac = self._cellularity * pop_frac
    vaf = tumor_cells_frac / 2
    ref_reads = int((1 - vaf) * total_reads)
    return ref_reads

  def _calc_total_reads(self, pop_frac, locus_start, locus_end, new_cn):
    # Proportion of all cells carrying CNV.
    p = self._cellularity * pop_frac
    if new_cn == 2:
      # If no net change in copy number -- e.g., because (major, minor) went
      # from (1, 1) to (2, 0) -- force the delta_cn to be 1.
      delta_cn = 1.
      no_net_change = True
    else:
      delta_cn = float(new_cn - 2)
      no_net_change = False

    region_length = locus_end - locus_start + 1
    fn = (self._read_depth * region_length) / self._read_length

    # This is a hack to prevent division by zero (when delta_cn = -2). Its
    # effect will be to make d large.
    if p == 1.0:
      p = 0.999

    d = (delta_cn**2 / 4) * (fn * p * (2 - p)) / (1 + (delta_cn  * p) / 2)

    if no_net_change:
      # If no net change in CN occurred, the estimate was just based on BAFs,
      # meaning we have lower confidence in it. Indicate this lack of
      # confidence via d by multiplying it by (read length / distance between
      # common SNPs), with the "distance between common SNPs" taken to be 1000 bp.
      d *= (self._read_length / 1000.)

    # Cap at 1e6 * read_depth.
    return int(round(min(d, self._max_reads())))

  def _format_overlapping_variants(self, variants, maj_cn, min_cn):
      variants = [(ssm_id, str(min_cn), str(maj_cn)) for ssm_id in variants]
      return variants

  def _format_cnvs(self, cnvs, variants):
    log('Estimated read depth: %s' % self._read_depth)

    for chrom, chrom_cnvs in cnvs.items():
      for cnv in chrom_cnvs:
        overlapping_variants = self._find_overlapping_variants(chrom, cnv, variants)
        total_reads = self._calc_total_reads(
          cnv['frac'],
          cnv['start'],
          cnv['end'],
          cnv['nmaj'] + cnv['nmin'],
        )
        yield {
          'chrom': chrom,
          'start': cnv['start'],
          'end': cnv['end'],
          'nmaj': cnv['nmaj'],
          'nmin': cnv['nmin'],
          'frac': cnv['frac'],
          'ref_reads': self._calc_ref_reads(cnv['frac'], total_reads),
          'total_reads': total_reads,
          'overlapping_variants': self._format_overlapping_variants(overlapping_variants, cnv['nmaj'], cnv['nmin']),
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
  def format_and_merge_cnvs(self, cnvs, variants):
    formatted = list(self._format_cnvs(cnvs, variants))
    formatted.sort(key = lambda f: f['frac'])
    if len(formatted) == 0:
      return []

    merged, formatted = formatted[:1], formatted[1:]
    merged[0]['cnv_id'] = 'c0'
    counter = 1

    for current in formatted:
      last = merged[-1]

      # Only merge CNVs if they're clonal. If they're subclonal, leave them
      # free to move around the tree.
      if current['frac'] == last['frac'] == 1.0:
        # Merge the CNVs.
        log('Merging %s_%s and %s_%s' % (current['chrom'], current['start'], last['chrom'], last['start']))
        last['total_reads'] = current['total_reads'] + last['total_reads']
        last['ref_reads'] = self._calc_ref_reads(last['frac'], last['total_reads'])
        self._merge_variants(last, current)
      else:
        # Do not merge the CNVs.
        current['cnv_id'] = 'c%s' % counter
        merged.append(current)
        counter += 1

    for cnv in merged:
      cnv['ref_reads'] = int(round(cnv['ref_reads'] * self._cnv_confidence))
      cnv['total_reads'] = int(round(cnv['total_reads'] * self._cnv_confidence))

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

  def format_variants(self, variant_list, error_rate):
    variant_list.sort(key = lambda v: variant_key(v[0]))

    for variant, ref_reads, total_reads in variant_list:
      ssm_id = 's%s' % self._counter
      variant_name = '%s_%s' % (variant.CHROM, variant.POS)

      # TODO: switch back to using calc_ref_freq() when we no longer want mu_r
      # and mu_v fixed.
      # This is mu_r in PhyloWGS.
      expected_ref_freq = 1 - error_rate
      if variant.CHROM in ('x', 'y', 'm'):
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
        'ref_reads': ref_reads,
        'total_reads': total_reads,
        'expected_ref_freq': expected_ref_freq,
        'expected_var_freq': expected_var_freq,
      }
      self._counter += 1

def restricted_float(x):
  x = float(x)
  if x < 0.0 or x > 1.0:
    raise argparse.ArgumentTypeError('%r not in range [0.0, 1.0]' % x)
  return x

def variant_key(var):
  chrom = var.CHROM
  if chrom == 'x':
    chrom = 100
  elif chrom == 'y':
    chrom = 101
  else:
    chrom = int(chrom)
  return (chrom, var.POS)

class VariantAndCnvGroup(object):
  def __init__(self):
    self._cn_regions = None

  def add_variants(self, variants_and_reads):
    self._variants_and_reads = variants_and_reads
    # Estimate read depth before any filtering of variants is performed, in
    # case no SSMs remain afterward.
    self._estimated_read_depth = self._estimate_read_depth()

  def add_cnvs(self, cn_regions):
    self._cn_regions = cn_regions

  def has_cnvs(self):
    return self._cn_regions is not None

  def _filter_variants_outside_regions(self, variants_and_reads, regions):
    filtered = []

    for variant, ref_reads, total_reads in self._variants_and_reads:
      for region in regions[variant.CHROM]:
        if region['start'] <= variant.POS <= region['end']:
          filtered.append((variant, ref_reads, total_reads))
          break

    return filtered

  def _is_region_normal_cn(self, region):
    return region['nmaj'] == region['nmin'] == 1

  def _print_variant_differences(self, before, after, before_label, after_label):
    before = set([var for (var, _, _) in before])
    after = set([var for (var, _, _) in after])
    log('%s=%s %s=%s delta=%s' % (before_label, len(before), after_label, len(after), len(before) - len(after)))


    assert after.issubset(before)
    removed = list(before - after)
    removed.sort(key = variant_key)

    for var in removed:
      var_name = '%s_%s' % (var.CHROM, var.POS)
      for region in self._cn_regions[var.CHROM]:
        if region['start'] <= var.POS <= region['end']:
          region_type = (self._is_region_normal_cn(region) and 'normal') or 'abnormal'
          log('%s\t[in %s-CN region chr%s(%s, %s)]' % (var_name, region_type, var.CHROM, region['start'], region['end']))
          break
      else:
        log('%s\t[outside all regions]' % var_name)

  def retain_only_variants_in_normal_cn_regions(self):
    if not self.has_cnvs():
      raise Exception('CN regions not yet provided')

    normal_cn = defaultdict(list)

    for chrom, regions in self._cn_regions.items():
      for region in regions:
        if self._is_region_normal_cn(region) and region['frac'] == 1:
          normal_cn[chrom].append(region)

    filtered = self._filter_variants_outside_regions(self._variants_and_reads, normal_cn)
    self._print_variant_differences(self._variants_and_reads, filtered, 'all_variants', 'only_normal_cn')
    self._variants_and_reads = filtered

  def _filter_multiple_abnormal_cn_regions(self, regions):
    good_regions = defaultdict(list)
    for chrom, reg in regions.items():
      idx = 0
      while idx < len(reg):
        region = reg[idx]

        # Accept clonal regions unconditonally, whether normal or abnormal CN.
        if region['frac'] == 1.0:
          good_regions[chrom].append(region)
          idx += 1

        else:
          regions_at_same_coords = [region]

          i = idx + 1
          while i < len(reg) and reg[i]['start'] == region['start']:
            # Battenerg either has entire regions at same coords, or they have
            # no overlap whatsoever. Thus, this assertion maintains sanity for
            # Battenberg, but may fail for other CN callers.
            assert reg[i]['end'] == region['end']
            regions_at_same_coords.append(reg[i])
            i += 1

          # This is assumption again specific to Battenberg.
          assert len(regions_at_same_coords) == 2
          abnormal_regions = [r for r in regions_at_same_coords if not self._is_region_normal_cn(r)]
          # In Battenberg, either one region is normal and the other abnormal, or both are abnormal.
          assert 1 <= len(abnormal_regions) <= 2
          # Ignore normal region and add only the single abnormal one. We do
          # this so PWGS can recalculate the frequencies based on the
          # major/minor CN of the region, according to the a & d values we will
          # assign to the region. (Does this make sense?)
          if len(abnormal_regions) == 1:
            good_regions[chrom].append(abnormal_regions[0])
          else:
            # Ignore CNV regions with multiple abnormal CN states.
            pass
          idx += len(regions_at_same_coords)

    return good_regions


  def exclude_variants_in_subclonal_cnvs(self):
    # Five possible placements for variant in Battenberg according to CN records:
    # 1 record:
    #   That record has normal CN: include
    #   That record has abnormal CN: include
    # 2 records:
    #   One record is normal CN, one record is abnormal CN: include
    #   Both records are abnormal CN: exclude (as we don't know what order the CN events occurred in)
    # Variant isn't listed in *any* region: exclude (as we suspect Battenberg didn't know what to do with the region)
    if not self.has_cnvs():
      raise Exception('CN regions not yet provided')

    good_regions = self._filter_multiple_abnormal_cn_regions(self._cn_regions)
    filtered = self._filter_variants_outside_regions(self._variants_and_reads, good_regions)
    self._print_variant_differences(self._variants_and_reads, filtered, 'all_variants', 'outside_subclonal_cn')
    self._variants_and_reads = filtered

  def subsample_variants(self, sample_size):
    random.shuffle(self._variants_and_reads)
    self._variants_and_reads = self._variants_and_reads[:sample_size]

  def write_variants(self, outfn, error_rate):
    formatter = VariantFormatter()
    self._variants = list(formatter.format_variants(self._variants_and_reads, error_rate))

    with open(outfn, 'w') as outf:
      print('\t'.join(('id', 'gene', 'a', 'd', 'mu_r', 'mu_v')), file=outf)
      for variant in self._variants:
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
    if len(self._variants_and_reads) == 0:
      return 50
    for variant, ref_reads, total_reads in self._variants_and_reads:
      read_sum += total_reads
    return float(read_sum) / len(self._variants_and_reads)

  def write_cnvs(self, outfn, cnv_confidence, cellularity, read_length):
    abnormal_regions = {}
    filtered_regions = self._filter_multiple_abnormal_cn_regions(self._cn_regions)
    for chrom, regions in filtered_regions.items():
      abnormal_regions[chrom] = [r for r in regions if not self._is_region_normal_cn(r)]

    with open(outfn, 'w') as outf:
      print('\t'.join(('cnv', 'a', 'd', 'ssms')), file=outf)
      formatter = CnvFormatter(cnv_confidence, cellularity, self._estimated_read_depth, read_length)
      for cnv in formatter.format_and_merge_cnvs(abnormal_regions, self._variants):
        overlapping = [','.join(o) for o in cnv['overlapping_variants']]
        vals = (
          cnv['cnv_id'],
          str(cnv['ref_reads']),
          str(cnv['total_reads']),
          ';'.join(overlapping),
        )
        print('\t'.join(vals), file=outf)

def log(msg):
  if log.verbose:
    print(msg)
log.verbose = False

def main():
  parser = argparse.ArgumentParser(
    description='Create ssm_dat.txt and cnv_data.txt input files for PhyloWGS from VCF and CNV data.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('-e', '--error-rate', dest='error_rate', type=restricted_float, default=0.001,
    help='Expected error rate of sequencing platform')
  parser.add_argument('-s', '--sample-size', dest='sample_size', type=int,
    help='Subsample SSMs to reduce PhyloWGS runtime')
  parser.add_argument('-b', '--battenberg', dest='battenberg',
    help='Path to Battenberg-formatted list of CNVs')
  parser.add_argument('--only-normal-cn', dest='only_normal_cn', action='store_true', default=False,
      help='Only output variants lying in normal CN regions. Do not output CNV data directly.')
  parser.add_argument('--output-cnvs', dest='output_cnvs', default='cnv_data.txt',
    help='Output destination for CNVs')
  parser.add_argument('--output-variants', dest='output_variants', default='ssm_data.txt',
    help='Output destination for variants')
  parser.add_argument('-c', '--cellularity', dest='cellularity', type=restricted_float, default=1.0,
    help='Fraction of sample that is cancerous rather than somatic. Used only for estimating CNV confidence -- if no CNVs, need not specify argument.')
  parser.add_argument('-v', '--variant-type', dest='input_type', required=True, choices=('sanger', 'mutect_pcawg', 'mutect_smchet','muse','dkfz','strelka', 'vardict'),
    help='Type of VCF file')
  parser.add_argument('--tumor-sample', dest='tumor_sample',
    help='Name of the tumor sample in the input VCF file. Defaults to last sample if not specified.')
  parser.add_argument('--cnv-confidence', dest='cnv_confidence', type=restricted_float, default=1.0,
    help='Confidence in CNVs. Set to < 1 to scale "d" values used in CNV output file')
  parser.add_argument('--read-length', dest='read_length', type=int, default=100,
    help='Approximate length of reads. Used to calculate confidence in CNV frequencies')
  parser.add_argument('--verbose', dest='verbose', action='store_true')
  parser.add_argument('--muse-tier', dest='muse_tier', type=int, default=0,
    help='Maximum MuSE tier to include')
  parser.add_argument('vcf_file')
  args = parser.parse_args()

  log.verbose = args.verbose

  # Fix random seed to ensure same set of SSMs chosen when subsampling on each
  # invocation.
  random.seed(1)

  grouper = VariantAndCnvGroup()
  if args.input_type == 'sanger':
    variant_parser = SangerParser(args.vcf_file)
  elif args.input_type == 'mutect_pcawg':
    variant_parser = MutectPcawgParser(args.vcf_file, args.tumor_sample)
  elif args.input_type == 'mutect_smchet':
    variant_parser = MutectSmchetParser(args.vcf_file, args.tumor_sample)
  elif args.input_type == 'muse':
    variant_parser = MuseParser(args.vcf_file, args.muse_tier)
  elif args.input_type == 'dkfz':
    variant_parser = DKFZParser(args.vcf_file)
  elif args.input_type == 'strelka':
    variant_parser = StrelkaParser(args.vcf_file)
  elif args.input_type == 'vardict':
    variant_parser = VarDictParser(args.vcf_file, args.tumor_sample)
  variants_and_reads = variant_parser.list_variants()
  grouper.add_variants(variants_and_reads)

  if args.battenberg:
    cnv_parser = BattenbergParser(args.battenberg)
    cn_regions = cnv_parser.parse()
    grouper.add_cnvs(cn_regions)

  if args.only_normal_cn:
    grouper.retain_only_variants_in_normal_cn_regions()
  elif grouper.has_cnvs():
    grouper.exclude_variants_in_subclonal_cnvs()

  if args.sample_size:
    grouper.subsample_variants(args.sample_size)

  grouper.write_variants(args.output_variants, args.error_rate)
  if not args.only_normal_cn and grouper.has_cnvs():
    grouper.write_cnvs(args.output_cnvs, args.cnv_confidence, args.cellularity, args.read_length)

if __name__ == '__main__':
  main()
