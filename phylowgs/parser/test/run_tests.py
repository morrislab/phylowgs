import subprocess
import tempfile
import filecmp
import os
import shutil

def compare(vcf_format, input_vcf, good_ssm_output):
  output_dir = tempfile.mkdtemp()
  ssm_output_path = os.path.join(output_dir, 'output.ssm')
  cnv_output_path = os.path.join(output_dir, 'output.cnv')
  cmd =[
    'python2',
    '../create_phylowgs_inputs.py',
    '--regions', 'all',
    '--output-variants', ssm_output_path,
    '--output-cnvs', cnv_output_path,
    '--vcf-type', 'S1=%s' % vcf_format,
    'S1=%s' % input_vcf,
  ]
  subprocess.call(cmd)

  variant_output_matches = filecmp.cmp(ssm_output_path, good_ssm_output)
  cnv_output_matches = os.path.getsize(cnv_output_path) == 0
  shutil.rmtree(output_dir)
  if variant_output_matches and cnv_output_matches:
    print('%s passed (%s)' % (vcf_format, input_vcf))
  else:
    raise Exception('Outputs do not match for %s' % vcf_format)

def main():
  #for vcf_format in ('dkfz', 'muse', 'mutect_pcawg', 'mutect_smchet', 'sanger', 'vardict'):
  for vcf_format in ('vardict', 'strelka', 'mutect_tcga'):
    compare(
      vcf_format,
      os.path.join('inputs', vcf_format, vcf_format + '.vcf'),
      os.path.join('outputs', vcf_format, vcf_format + '.ssm')
    )
  return

  # test for ID parsing
  vcf_format = 'mutect_tcga'
  compare(
    vcf_format,
    os.path.join('inputs', vcf_format, vcf_format + '_with_ids.vcf'),
    os.path.join('outputs', vcf_format, vcf_format + '_with_ids.ssm')
  )

main()
