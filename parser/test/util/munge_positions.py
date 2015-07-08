from collections import defaultdict
import random
import sys

def chrom_key(chrom):
  chrom = chrom.lower()
  if chrom.isdigit():
    return int(chrom)
  elif chrom == 'x':
    return 101
  elif chrom == 'y':
    return 102
  else:
    return 999

def main():
  variants = defaultdict(list)

  for line in sys.stdin:
    line = line.strip()
    if line.startswith('#'):
      print(line)
      continue

    fields = line.split('\t')
    chrom = fields[0]
    variants[chrom].append(fields)

  for chrom in sorted(variants.keys(), key = chrom_key):
    chrom_vars = variants[chrom]
    random.shuffle(chrom_vars)
    last_pos = 1

    # Make our chromosome span approximately 200 Mb. Given uniform distribution
    # of randomly drawn ints, half will result in smaller step size than needed
    # to fill 200 Mb, while half will be bigger than necessary. These forces
    # should approximately balance.
    max_step_size = int(2 * (200e6 / float(len(chrom_vars))))

    for var in chrom_vars:
      new_pos = last_pos + random.randint(1, max_step_size)
      last_pos = new_pos
      var[1] = str(new_pos)
      print('\t'.join(var))

main()
