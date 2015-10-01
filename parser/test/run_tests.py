import subprocess
import tempfile
import filecmp
import os
import shutil

def compare(vcf_format, input_vcf, good_ssm_output):
  output_dir = tempfile.mkdtemp()
  ssm_output_path = os.path.join(output_dir, 'output.ssm')
  cnv_output_path = os.path.join(output_dir, 'output.cnv')
  subprocess.call([
    'python2',
    '../create_phylowgs_inputs.py',
    '-v', vcf_format,
    '--output-variants', ssm_output_path,
    '--output-cnvs', cnv_output_path,
    input_vcf
  ])

  variant_output_matches = filecmp.cmp(ssm_output_path, good_ssm_output)
  cnv_output_matches = os.path.getsize(cnv_output_path) == 0
  shutil.rmtree(output_dir)
  if variant_output_matches and cnv_output_matches:
    print('%s passed' % vcf_format)
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

main()
