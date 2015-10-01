import numpy as np
from collections import defaultdict

class ResultMunger(object):
  def __init__(self, tree_summaries, mutlist, mutass, min_ssms):
    self._tree_summaries = tree_summaries
    self._mutlist = mutlist
    self._mutass = mutass
    self._min_ssms = min_ssms

    self._subclone_idx_map = defaultdict(dict)

  def _convert_keys_to_ints(self, dic):
    keys = dic.keys()
    for key in dic.keys():
      dic[int(key)] = dic[key]
      del dic[key]

  def munge(self):
    for tree_idx, tree_features in self._tree_summaries.items():
      small_nodes = self._remove_small_nodes(tree_idx, tree_features['populations'])
      self._remove_nodes_from_tree_structure(small_nodes, tree_features['structure'])
      self._renumber_nodes(tree_idx)

      self._reassign_muts(tree_idx)
      # Correct mutation counts.
      for sidx, subclone in tree_features['populations'].items():
        # Note that only mutass entries for subclones with (> 0 SSMs or > 0
        # CNVs) will exist. Thus, no mutass entry will exist for node 0, as it
        # never has SSMs or CNVs.
        if sidx == 0:
          continue
        for mut_type in ('ssms', 'cnvs'):
          subclone['num_%s' % mut_type] = len(self._mutass[tree_idx][sidx][mut_type])

    return (self._tree_summaries, self._mutass)

  def _renumber_nodes(self, tree_idx):
    subclone_idxs = sorted(self._tree_summaries[tree_idx]['populations'].keys())

    num_removed = 0
    # We may have removed populations beyond max(subclone_idxs), but as these
    # occurred *after* the highest-indexed of the remaining populations, so
    # renumbering is necessary for them.
    for subclone_idx in range(1, max(subclone_idxs) + 1):
      if subclone_idx not in self._tree_summaries[tree_idx]['populations']:
        num_removed += 1
      elif num_removed > 0:
        self._subclone_idx_map[tree_idx][subclone_idx] = subclone_idx - num_removed

    # By proceeding in sorted order, we guarantee we're not overwriting a
    # single element twice, which would give the wrong values.
    for subclone_idx in subclone_idxs:
      if subclone_idx not in self._subclone_idx_map[tree_idx]:
        # Node was removed, so do nothing.
        continue
      new_idx = self._subclone_idx_map[tree_idx][subclone_idx]

      self._tree_summaries[tree_idx]['populations'][new_idx] = self._tree_summaries[tree_idx]['populations'][subclone_idx]
      del self._tree_summaries[tree_idx]['populations'][subclone_idx]

      if subclone_idx in self._tree_summaries[tree_idx]['structure']:
        self._tree_summaries[tree_idx]['structure'][new_idx] = self._tree_summaries[tree_idx]['structure'][subclone_idx]
        del self._tree_summaries[tree_idx]['structure'][subclone_idx]

    # We must also renumber children in the structure -- just renumbering
    # parents isn't enough.
    for subclone_idx, children in self._tree_summaries[tree_idx]['structure'].items():
      self._tree_summaries[tree_idx]['structure'][subclone_idx] = [
        self._subclone_idx_map[tree_idx][c]
        if c in self._subclone_idx_map[tree_idx]
        else c
        for c in children
      ]

  def _remove_small_nodes(self, tree_idx, populations):
    small_nodes = set()

    subclone_idxs = sorted(populations.keys())
    last_phi = None
    last_idx = None

    if self._min_ssms >= 1:
      # This is a count of SSMs, so use it without adjustment (but ensure it's an int).
      min_ssms = int(self._min_ssms)
    else:
      # This is a fraction of total SSMs.
      num_ssms = len(self._mutlist['ssms'])
      min_ssms = int(round(self._min_ssms * num_ssms))

    for subclone_idx in subclone_idxs:
      for p, children in self._tree_summaries[tree_idx]['structure'].items():
        if subclone_idx in children:
          parent = p
          break
      subclone = populations[subclone_idx]
      # Ensure this node's phi is <= the phi of its preceding sibling node, if any exists.
      if subclone_idx > 0 and last_idx in self._tree_summaries[tree_idx]['structure'][parent]:
        assert np.mean(subclone['cellular_prevalence']) <= last_phi
      last_phi = np.mean(subclone['cellular_prevalence'])
      last_idx = subclone_idx

      if subclone_idx == 0 or subclone['num_ssms'] >= min_ssms:
        continue
      small_nodes.add(subclone_idx)

    for node_idx in small_nodes:
      del populations[node_idx]
      # Mark node as removed.
      self._subclone_idx_map[tree_idx][node_idx] = None

    return small_nodes

  def _remove_nodes_from_tree_structure(self, removed, tree_structure):
    def _find_parent(struct, idx):
      for parent, children in struct.items():
        if idx in children:
          return parent
      raise Exception('Could not find parent of %s in %s' % (idx, struct))

    removed = set(removed)
    for rem in removed:
      parent = _find_parent(tree_structure, rem)
      # Remove node from parent
      tree_structure[parent] = [c for c in tree_structure[parent] if c != rem]
      # Assign removed node's children to their grandparent
      if rem in tree_structure:
        tree_structure[parent] += tree_structure[rem]
        del tree_structure[rem]
      # Sort, since order may not be preserved.
      tree_structure[parent].sort()
      # If no children remain after deletion, remove child list from tree.
      if len(tree_structure[parent]) == 0:
        del tree_structure[parent]

  def _move_muts_to_best_node(self, muts, mutass, populations):
    for mut_type in ('ssms', 'cnvs'):
      for mut_id in muts[mut_type]:
        mut_stats = self._mutlist[mut_type][mut_id]
        ref_reads = np.mean(mut_stats['ref_reads'])
        total_reads = np.mean(mut_stats['total_reads'])
        # Note this doesn't take into account CNVs that may skew the
        # relationship between VAF and phi -- we assume that there is one
        # maternal and one paternal copy, and that only one of these is
        # mutated.
        implied_phi = 2 * (total_reads - ref_reads) / total_reads
        implied_phi = min(implied_phi, 1.0)

        lowest_phi_delta = 1
        best_node = None
        for pidx, pop in populations.items():
          phi_delta = abs(np.mean(pop['cellular_prevalence']) - implied_phi)
          # Don't allow assignments to the non-cancerous root node.
          if phi_delta < lowest_phi_delta and pidx != 0:
            lowest_phi_delta = phi_delta
            best_node = pidx

        mutass[best_node][mut_type].append(mut_id)

  def _reassign_muts(self, tree_idx):
    deleted_muts = []
    mutass = self._mutass[tree_idx]

    for sidx in sorted(self._subclone_idx_map[tree_idx].keys()):
      new_idx = self._subclone_idx_map[tree_idx][sidx]
      # This ensures we're not improperly overwriting assignments.
      assert new_idx < sidx

      if new_idx is None: # Node was removed.
        # This will comprise both SSMs and CNVs.
        deleted_muts.append(mutass[sidx])
      else:
        mutass[new_idx] = mutass[sidx]
      del mutass[sidx]

    for dm in deleted_muts:
      self._move_muts_to_best_node(dm, mutass, self._tree_summaries[tree_idx]['populations'])
