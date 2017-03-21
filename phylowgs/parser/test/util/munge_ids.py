import sys
import uuid

def main():
  new_id = uuid.uuid4()

  for line in sys.stdin:
    line = line.strip()
    if line.startswith('#'):
      print(line)
      continue

    fields = line.split('\t')
    fields[2] = str(new_id)
    print('\t'.join(fields))

main()
