import os
import glob
import json
from collections import defaultdict

def main():
  datasets = defaultdict(list)

  base_dir = 'data'
  for run_name in os.listdir(base_dir):
    for dataset_path in glob.glob(os.path.join(base_dir, run_name, '*.summ.json.gz')):
      dataset_name = dataset_path.split('/')[-1].split('.')[0]
      datasets[run_name].append({
        'path': dataset_path,
        'name': dataset_name,
      })

  print(json.dumps(datasets))

main()
