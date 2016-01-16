import gzip
import json
import zipfile

class ResultLoader(object):
  def __init__(self, tree_summary_fn, mutation_list_fn, mutation_assignment_fn):
    self._tree_summary_fn = tree_summary_fn
    self._mutation_list_fn = mutation_list_fn
    self._mutation_assignment_fn = mutation_assignment_fn

    self.mutlist = None
    self.tree_summary = None
    self.dataset_name = None

    self._load_tree_data()

  def _convert_keys_to_ints(self, dic):
    keys = dic.keys()
    for key in dic.keys():
      dic[int(key)] = dic[key]
      del dic[key]

  def _load_tree_data(self):
    with gzip.GzipFile(self._tree_summary_fn) as treesummf:
      tree_json = json.load(treesummf)
      self.dataset_name = tree_json['dataset_name']
      self.tree_summary = tree_json['trees']

    self._convert_keys_to_ints(self.tree_summary)
    for tree_idx, tree_features in self.tree_summary.items():
      self._convert_keys_to_ints(tree_features['populations'])
      self._convert_keys_to_ints(tree_features['structure'])

    with gzip.GzipFile(self._mutation_list_fn) as mutlistf:
      self.mutlist = json.load(mutlistf)
    self.num_ssms = len(self.mutlist['ssms'])

  def _load_assignments(self, mutf, tree_idx):
    mutass = json.loads(mutf.read('%s.json' % tree_idx))
    mutass = mutass['mut_assignments']
    self._convert_keys_to_ints(mutass)
    return mutass

  def load_mut_assignments(self, tree_idx):
    with zipfile.ZipFile(self._mutation_assignment_fn) as mutf:
      return self._load_assignments(mutf, tree_idx)

  def load_all_mut_assignments(self):
    with zipfile.ZipFile(self._mutation_assignment_fn) as mutf:
      tree_indices = [int(i.filename.split('.')[0]) for i in mutf.infolist()]
      tree_indices.sort()
      for tree_idx in tree_indices:
        yield (tree_idx, self._load_assignments(mutf, tree_idx))
