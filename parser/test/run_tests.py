import subprocess
import tempfile
import filecmp
import os
import shutil

def compare_outputs(comparison, parser_params=None):
  if parser_params is None:
    parser_params = []
  if 'cnv_path' not in comparison:
    comparison['cnv_path'] = []

  output_dir = tempfile.mkdtemp()
  ssm_output_path = os.path.join(output_dir, 'output.ssm')
  cnv_output_path = os.path.join(output_dir, 'output.cnv')
  cmd = [
    'python2',
    '../create_phylowgs_inputs.py',
    '--output-variants', ssm_output_path,
    '--output-cnvs', cnv_output_path,
  ]
  cmd += parser_params
  for idx, vcf_format in enumerate(comparison['vcf_format']):
    samp_name = 'S%s' % (idx + 1)
    cmd += ['--vcf-type', '%s=%s' % (samp_name, vcf_format)]
  for idx, cnv_path in enumerate(comparison['cnv_path']):
    samp_name = 'S%s' % (idx + 1)
    cmd += ['--cnvs', '%s=%s' % (samp_name, cnv_path)]
  for idx, vcf_path in enumerate(comparison['vcf_path']):
    samp_name = 'S%s' % (idx + 1)
    cmd.append('%s=%s' % (samp_name, vcf_path))

  subprocess.call(cmd)

  files = (
    (ssm_output_path, comparison['ssm_output_good']),
    (cnv_output_path, comparison['cnv_output_good']),
  )
  all_match = True
  for inf, outf in files:
    all_match = all_match and filecmp.cmp(inf, outf)
  shutil.rmtree(output_dir)
  if all_match:
    print('%s passed' % comparison['name'])
  else:
    raise Exception('Outputs do not match for %s' % comparison['name'])

def test_vcf_formats():
  for vcf_format in ('vardict', 'strelka', 'mutect_tcga'):
    comparison = {
      'name':            'vcf_format_%s' % vcf_format,
      'vcf_format':      [vcf_format],
      'vcf_path':        [os.path.join('inputs',  vcf_format, vcf_format + '.vcf')],
      'ssm_output_good': os.path.join('outputs', vcf_format, vcf_format + '.ssm'),
      'cnv_output_good': os.path.join('outputs', vcf_format, vcf_format + '.cnv'),
    }
    compare_outputs(comparison, ['--regions', 'all'])

def test_id_parsing():
  vcf_format = 'mutect_tcga'
  comparison = {
    'name':            'id_parsing',
    'vcf_format':      [vcf_format],
    'vcf_path':        [os.path.join('inputs',  vcf_format, vcf_format + '_with_ids.vcf')],
    'ssm_output_good': os.path.join('outputs', vcf_format, vcf_format + '_with_ids.ssm'),
    'cnv_output_good': os.path.join('outputs', vcf_format, vcf_format + '_with_ids.cnv'),
  }
  compare_outputs(comparison, ['--regions', 'all'])

def test_singlesamp_cnv():
  comparison = {
    'name':            'singlesamp_cnv',
    'vcf_format':      ['pcawg_consensus'],
    'vcf_path':        [os.path.join('inputs', 'singlesamp_cnv/S1.vcf.gz')],
    'cnv_path':        [os.path.join('inputs', 'singlesamp_cnv/S1.cnv.txt')],
    'ssm_output_good': os.path.join('outputs', 'singlesamp_cnv', 'ssm_data.txt'),
    'cnv_output_good': os.path.join('outputs', 'singlesamp_cnv', 'cnv_data.txt'),
  }
  compare_outputs(comparison)

def test_multisamp_cnvs():
  comparison = {
    'name':            'multisamp_cnvs',
    'vcf_format':      2*['pcawg_consensus'],
    'vcf_path':        [os.path.join('inputs', 'multisamp_cnvs/S%s.vcf.gz') % sidx for sidx in (1, 2)],
    'cnv_path':        [os.path.join('inputs', 'multisamp_cnvs/S%s.cnv.txt') % sidx for sidx in (1, 2)],
    'ssm_output_good': os.path.join('outputs', 'multisamp_cnvs', 'ssm_data.txt'),
    'cnv_output_good': os.path.join('outputs', 'multisamp_cnvs', 'cnv_data.txt'),
  }
  compare_outputs(comparison)

def main():
  test_vcf_formats()
  test_singlesamp_cnv()
  test_multisamp_cnvs()

  # This test is broken right now. D'oh.
  #test_id_parsing()

main()
