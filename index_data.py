import os
import glob
import json
from collections import defaultdict

def main():
  datasets = defaultdict(list)

  base_dir = 'data'
  for run_name in os.listdir(base_dir):
    for summary_path in glob.glob(os.path.join(base_dir, run_name, '*.summ.json')):
      dataset_name = summary_path.split('/')[-1].split('.')[0]
      datasets[run_name].append({
        'summary_path': summary_path,
        'muts_path': os.path.join(base_dir, run_name, dataset_name + '.muts.json'),
        'name': dataset_name,
      })

  print(json.dumps(datasets))

main()
