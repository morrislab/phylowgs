import json
import gzip
import zipfile

class JsonWriter(object):
  def __init__(self, dataset_name):
    self._dataset_name = dataset_name

  def write_mutlist(self, mutlist, mutlist_outfn):
    with gzip.GzipFile(mutlist_outfn, 'w') as mutf:
      mutlist['dataset_name'] = self._dataset_name
      json.dump(mutlist, mutf)

  def write_summaries(self, summaries, params, summaries_outfn):
    to_dump = {
      'dataset_name': self._dataset_name,
      'params': params,
      'trees': summaries,
    }
    with gzip.GzipFile(summaries_outfn, 'w') as summf:
      json.dump(to_dump, summf)

  def write_mutass(self, mutass, mutass_outfn):
    with zipfile.ZipFile(mutass_outfn, 'w', compression=zipfile.ZIP_DEFLATED) as muts_file:
      for tree_idx, tree_mutass in mutass.items():
        to_dump = {
          'mut_assignments': tree_mutass,
          'dataset_name': self._dataset_name
        }
        muts_file.writestr('%s.json' % tree_idx, json.dumps(to_dump))

