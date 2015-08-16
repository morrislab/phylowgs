import argparse
import util2
import util
import tssb
import numpy as np
import data
import json
import os

def post_assignments(name, id, a, d, mu_r, mu_v, tree_file, cnv_ids=None, copies=None):
  reader = util2.TreeReader(tree_file)
  assignments = []

  for idx, llh, tree in reader.load_trees_and_metadata(remove_empty_vertices = True):
    data_ob = data.Datum(name, id, a, d, mu_r, mu_v)
    data_ob.tssb = tree
    if cnv_ids:
      cnv_obs = [x for x in tree.data if x.name in cnv_ids]
      for i,c in enumerate(cnv_obs):
        data_ob.cnv.append((c,copies[i][0],copies[i][1]))
      
    nodes = tree.get_nodes()
    tree.data.append(data_ob)
    tree.assignments.append(nodes[1])
    nodes[1].add_datum(len(tree.data)-1)
    probs = np.array([compute_llh(n,tree,data_ob) for n in nodes[1:]])
    probs = np.exp(probs - util.logsumexp(probs))
    node_index = sum( np.random.rand() > np.cumsum(probs) ) + 1
    mapping = construct_index_map(tree,nodes)
    assignments.append(mapping[node_index])

  return assignments

def construct_index_map(tree, nodes):
  idx = [0]
  mapping = {}
  def decend(node):
    mapping[nodes.index(node)] = idx[0]
    children = sorted(node.children(), key = lambda x: np.mean(x.params), reverse=True)
    for child in children:
      idx[0] += 1
      decend(child)

  decend(tree.root['node'])
  return mapping

def compute_llh(node,tree,datum):
  n = len(tree.data)-1
  tree.assignments[-1].remove_datum(n)
  tree.assignments[-1] = node
  node.add_datum(n)
  return node.logprob(tree.data[n:n+1])

def read_ssms(ssm_file):
  with open(ssm_file) as ssmf:
    header = next(ssmf)
    for line in ssmf:
      fields = line.strip().split('\t')
      yield fields

def parse_cnvs(cnvs_file):
  cnvs = {}

  with open(cnvs_file) as cnvf:
    header = next(cnvf)
    for line in cnvf:
      fields = line.strip().split()
      if len(fields) == 4:
        cnv_id = fields[0]
        ssms = fields[3]
        overlapping_ssms = [(ssm_id, int(cn1), int(cn2)) for ssm_id, cn1, cn2 in [s.split(',') for s in ssms.split(';')]]
        cnvs[cnv_id] = overlapping_ssms
      else: # No overlapping SSMs provided.
        pass

  return cnvs

def find_overlapping_cnvs(ssm_id, cnvs):
  overlapping_cnvs = []

  for cnv_id, overlapping_ssms in cnvs.items():
    for sid, cn1, cn2 in overlapping_ssms:
      if sid != ssm_id:
        continue
      overlapping_cnvs.append((cnv_id, cn1, cn2))

  return overlapping_cnvs

def main():
  parser = argparse.ArgumentParser(
    description='Assign mutations to existing trees',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--cnvs', dest='cnv_file',
    help='Path to CNV list')
  parser.add_argument('ssm_file',
      help='SSM file whose mutations you wish to reassign')
  parser.add_argument('trees_file',
    help='Path to sampled trees')
  parser.add_argument('ssm_ids', nargs='+',
      help='ID of SSM in ssm_file to reassign')
  args = parser.parse_args()

  if args.cnv_file:
    cnvs = parse_cnvs(args.cnv_file)
  else:
    cnvs = {}

  ssm_assignments = {}

  for ssm_id, ssm_name, a, d, mu_r, mu_v in read_ssms(args.ssm_file):
    if ssm_id not in args.ssm_ids:
      continue
    overlapping_cnvs = find_overlapping_cnvs(ssm_id, cnvs)
    if overlapping_cnvs:
      cnv_ids = [c[0] for c in overlapping_cnvs]
      cnv_copies = [c[1:3] for c in overlapping_cnvs]
    else:
      cnv_ids = None
      cnv_copies = None

    a = [int(i) for i in a.split(',')]
    d = [int(i) for i in d.split(',')]
    mu_r = float(mu_r)
    mu_v = float(mu_v)

    ssm_assignments[ssm_id] = post_assignments(ssm_name, ssm_id, a, d, mu_r, mu_v, args.trees_file, cnv_ids, cnv_copies)

  print(json.dumps({
    'assignments': ssm_assignments,
    'trees': os.path.realpath(args.trees_file)
  }))


if __name__ == '__main__':
  #post_assignments('s1313', 's1313', [25], [50], .999, .5, 'trees.cnv.zip')
  #post_assignments('s1313', 's1313', [45], [50], .999, .5, 'trees.zip',cnv_ids=['c0'],copies=[(0,2)])
  main()
