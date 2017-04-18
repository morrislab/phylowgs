from __future__ import print_function
import numpy as np
import itertools

class IndexCalculator(object):
  def __init__(self, tree_summ):
    self._tree_summ = tree_summ
    self._poprel = self._determine_pop_relations(tree_summ['structure'])

  def _determine_pop_relations(self, tree_struct):
    relations = {}
    all_verts = set()

    def _traverse_r(vertex, ancestors):
      all_verts.add(vertex)
      for anc in ancestors:
        relations[(anc, vertex)] = 'anc_desc'
        relations[(vertex, anc)] = 'desc_anc'
      if vertex in tree_struct:
        for child in tree_struct[vertex]:
          _traverse_r(child, ancestors + [vertex])

    root = 0
    _traverse_r(root, [])

    for vert1, vert2 in itertools.combinations(all_verts, 2):
      if (vert1, vert2) in relations:
        continue
      relations[(vert1, vert2)] = 'cousin'
      relations[(vert2, vert1)] = 'cousin'

    return relations

  def _calc_index(self, reltype):
    tree_pops = self._tree_summ['populations']
    totalssms = sum([P['num_ssms'] for P in tree_pops.values()])
    index = 0

    for (popidx1, popidx2), relation in self._poprel.items():
      if relation != reltype:
        continue
      nssms1, nssms2 = tree_pops[popidx1]['num_ssms'], tree_pops[popidx2]['num_ssms']
      index += nssms1 * nssms2

    # The maximum value of `index` will be `N(N - 1)`, as we exclude the
    # diagonal. If we ignore the diagonal, then if we were to explicitly
    # calculate the N*N matrices corresponding to each of the indices and then
    # sum the matrices, every entry in the sum would be 1.
    normidx = float(index) / (totalssms * (totalssms - 1))
    assert 0. <= normidx <= 1.
    return normidx

  def calc_linearity_index(self):
    lowertrisum = self._calc_index('anc_desc')
    # Note that the matrix isn't symmetric -- if A is an ancestor of B, then we
    # know that B is *not* an ancestor of A. The clustering and branching
    # matrices are, however, symmetric. Thus, for this, we should be
    # normalizing against (N choose 2) instead of (N permute 2); but since
    # _calc_index() normalizes against the latter, we multiply by two to
    # correct for this, since (N permute 2) = 2(N choose 2).
    linidx = 2 * lowertrisum
    assert 0. <= linidx <= 1.
    return linidx

  def calc_branching_index(self):
    return self._calc_index('cousin')

  def calc_clustering_index(self):
    # Technically, in the clustering matrix, the diagonal should be 1's
    # (since a mutation should be said to cluster with itself). But to keep the
    # value in the denominator against which we normalize the same across all
    # three indices, thereby ensuring that the sum of the three normalized
    # indices will be 1, we force the diagonal to be zero -- which is why we
    # add nssms(nssms - 1) rather than nssms^2 below.
    ccidx = 0
    totalssms = 0
    for pop in self._tree_summ['populations'].values():
      nssms = pop['num_ssms']
      ccidx += nssms * (nssms - 1)
      totalssms += nssms
    normccidx = float(ccidx) / (totalssms * (totalssms - 1))

    assert 0. <= normccidx <= 1.
    # It's enough just to subtract the other two indices from one to calculate
    # the branching index. However, calculating it independently and then
    # checking it against that result gives more confidnece it's correct.
    assert np.isclose(normccidx, 1 - self.calc_linearity_index() - self.calc_branching_index())
    return normccidx
