import numpy as np
from collections import defaultdict

import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import util2

class ResultGenerator(object):
  def generate(self, tree_file, include_ssm_names):
    reader = util2.TreeReader(tree_file)
    first_tree = next(reader.load_trees())
    reader.close()
    mutlist = self._list_mutations(first_tree, include_ssm_names)

    summaries = {}
    all_mutass = {}
    for idx, llh, pops, mutass, structure in self._summarize_all_pops(tree_file):
      summaries[idx] = {
        'llh': llh,
        'structure': structure,
        'populations': pops,
      }
      all_mutass[idx] = mutass

    return summaries, mutlist, all_mutass

  def _summarize_all_pops(self, tree_file):
    reader = util2.TreeReader(tree_file)
    for idx, llh, tree in reader.load_trees_and_metadata(remove_empty_vertices = True):
      yield (idx, llh) + self._summarize_pops(tree)
    reader.close()

  def _summarize_pops(self, tree):
    pops = {}
    structure = defaultdict(list)
    # Note that there will be an entry in mut_assignments for a given subclone
    # only if it has at least one SSM or CNV. This assumption holds true for all
    # PhyloWGS trees, so one can ascertain the number of cancerous populations
    # via len(mut_assignments).
    mut_assignments = defaultdict(lambda: {'cnvs': [], 'ssms': []})
    idx = [0]

    def _traverse_r(vertex, parent):
      mutations = vertex.get_data()
      # vertex.params represents phis (i.e., population freqs) associated with
      # each sample.
      cell_prev = list(vertex.params)
      current_idx = idx[0]

      num_ssms = 0
      num_cnvs = 0
      for mut in mutations:
        if mut.id.startswith('s'):
          mut_assignments[current_idx]['ssms'].append(mut.id)
          num_ssms += 1
        elif mut.id.startswith('c'):
          mut_assignments[current_idx]['cnvs'].append(mut.id)
          num_cnvs += 1
        else:
          raise Exception('Unknown mutation ID type: %s' % mut.id)

      # Preorder traversal is consistent with printo_latex.py, meaning index
      # values should correspond to same vertices.
      pops[current_idx] = {
        'cellular_prevalence': cell_prev,
        'num_ssms': num_ssms,
        'num_cnvs': num_cnvs,
      }

      # Visit children in order of decreasing phi.
      children = sorted(vertex.children(), key = lambda v: np.mean(v.params), reverse=True)
      for child in children:
        idx[0] += 1
        structure[current_idx].append(idx[0])
        _traverse_r(child, current_idx)

    _traverse_r(tree.root['node'], None)
    return (pops, mut_assignments, structure)

  def _list_mutations(self, tree, include_ssm_names):
    cnvs = {}
    ssms = {}
    ssms_in_cnvs = defaultdict(list)

    def _traverse(node):
      for mut in node['node'].get_data():
        if mut.id.startswith('s'):
          ssms[mut.id] = {
            'ref_reads': mut.a,
            'total_reads': mut.d,
            'expected_ref_in_ref': mut.mu_r,
            'expected_ref_in_variant': mut.mu_v
          }
          if include_ssm_names:
            ssms[mut.id]['name'] = mut.name

          for cnv, maternal_cn, paternal_cn in mut.cnv:
            ssms_in_cnvs[cnv.id].append({
              'ssm_id': mut.id,
              'maternal_cn': maternal_cn,
              'paternal_cn': paternal_cn,
            })
        elif mut.id.startswith('c'):
          cnvs[mut.id] = {
            'ref_reads': mut.a,
            'total_reads': mut.d
          }
        else:
          raise Exception('Unknown mutation type: %s' % mut.id)
      for child in node['children']:
        _traverse(child)
    _traverse(tree.root)

    for cnv_id, cnv in cnvs.items():
      cnv['ssms'] = ssms_in_cnvs[cnv_id]

    return {
      'ssms': ssms,
      'cnvs': cnvs,
    }
