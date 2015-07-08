import subprocess
import tempfile
import filecmp
import os
import shutil

def compare(vcf_format, input_vcf, good_ssm_output):
  output_dir = tempfile.mkdtemp()
  output_path = os.path.join(output_dir, 'output.ssm')
  subprocess.call([
    'python2',
    '../create_phylowgs_inputs.py',
    '-v', vcf_format,
    '--output-variants', output_path,
    input_vcf
  ])

  output_matches = filecmp.cmp(output_path, good_ssm_output)
  shutil.rmtree(output_dir)
  if output_matches:
    print('%s passed' % vcf_format)
  else:
    raise Exception('VCFs do not match for %s' % vcf_format)

def main():
  #for vcf_format in ('dkfz', 'muse', 'mutect_pcawg', 'mutect_smchet', 'sanger', 'vardict'):
  for vcf_format in ('vardict',):
    compare(
      vcf_format,
      os.path.join('inputs', vcf_format, vcf_format + '.vcf'),
      os.path.join('outputs', vcf_format, vcf_format + '.ssm')
    )

main()
