#!/usr/bin/env python2
import argparse
import util2
import util
import tssb
import numpy as np
import data
import collections
import json
import os

def post_assignments(ssms, tree_file):
  epsilon = np.finfo(np.float64).eps
  reader = util2.TreeReader(tree_file)
  other_ssms = collections.defaultdict(list)
  assignments = collections.defaultdict(list)
  old_tree = None
  error = 0
  for idx, tree in reader.load_trees_and_burnin(remove_empty_vertices = True):
    tree.assignments.append(None)
    tree.data.append(None)
    for name, id, a, d, mu_r, mu_v, cnv_ids, copies in ssms:
      data_ob = data.Datum(name, id, a, d, mu_r, mu_v)
      data_ob.tssb = tree
      if cnv_ids:
        cnv_obs = [x for x in tree.data[:-1] if x.name in cnv_ids]
        for i,c in enumerate(cnv_obs):
          data_ob.cnv.append((c,copies[i][0],copies[i][1]))
      nodes = tree.get_nodes()
      tree.data[-1] = data_ob

      new_node = get_new_node(tree,other_ssms[name])

      tree.assignments[-1] = new_node
      n = len(tree.data)-1
      new_node.add_datum(n)
      llhmap = {}

      max_u = 1.0
      min_u = 0.0
      old_llh = tree.assignments[n].logprob(tree.data[n:n+1])
      llhmap[tree.assignments[n]] = old_llh
      llh_s = np.log(np.random.rand()) + old_llh

      while True:
        new_u = (max_u-min_u)*np.random.rand() + min_u
        new_node = find_node2(n,nodes,new_u)
        old_node = tree.assignments[n]                 
        old_node.remove_datum(n)
        new_node.add_datum(n)
        tree.assignments[n] = new_node
        if new_node in llhmap:
          new_llh = llhmap[new_node]
        else:
          new_llh = new_node.logprob(tree.data[n:n+1])
          llhmap[new_node] = new_llh
        if new_llh > llh_s:
          break
        elif abs(max_u-min_u) < epsilon:
          new_node.remove_datum(n)
          old_node.add_datum(n)
          tree.assignments[n] = old_node
          new_node = old_node
          #print >>sys.stderr, "Slice sampler shrank down.  Keep current state."
          break
        else:
          new_node.remove_datum(n)
          old_node.add_datum(n)
          tree.assignments[n] = old_node
          if nodes.index(old_node) > nodes.index(new_node):
            min_u = new_u
          else:
            max_u = new_u
    
      other_ssms[name] = list(new_node.data)
      other_ssms[name].remove(n)

      nodes = tree.get_nodes()
      mapping = construct_index_map(tree,nodes)
      assignments[id].append(mapping[nodes.index(new_node)])
  return assignments

def get_new_node(tree,other_ssms):
  if not other_ssms:
    return tree.get_nodes()[1]
  counts = collections.Counter([tree.assignments[i] for i in other_ssms])

  index = sum( np.random.rand() > np.cumsum([x/float(len(other_ssms)) for x in counts.values()]) )
  return counts.keys()[ sum( np.random.rand() > np.cumsum([x/float(len(other_ssms)) for x in counts.values()]) )]

def path_lt(path1, path2):
  if len(path1) == 0 and len(path2) == 0:
    return 0
  elif len(path1) == 0:
    return 1
  elif len(path2) == 0:
    return -1
  s1 = "".join(map(lambda i: "%03d" % (i), path1))
  s2 = "".join(map(lambda i: "%03d" % (i), path2))

  return cmp(s2, s1)

def find_node2(n,nodes,u):
  ndata = float(n)
  probs = np.array([n.num_local_data()/ndata for n in nodes])
  probs = u > np.cumsum(probs)
  index = sum(probs)
  return nodes[index]

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

  ssms = []
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
    ssms.append((ssm_name, ssm_id, a, d, mu_r, mu_v, cnv_ids, cnv_copies))

  ssm_assignments = post_assignments(ssms, args.trees_file)

  print(json.dumps({
    'assignments': ssm_assignments,
    'trees': os.path.realpath(args.trees_file)
  }))


if __name__ == '__main__':
  #post_assignments('s1313', 's1313', [25], [50], .999, .5, 'trees.cnv.zip')
  #post_assignments('s1313', 's1313', [45], [50], .999, .5, 'trees.zip',cnv_ids=['c0'],copies=[(0,2)])
  main()
