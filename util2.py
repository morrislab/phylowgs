import numpy
from numpy import *
import cPickle as pickle
import zipfile
import shutil

import scipy.stats as stat
from scipy.stats import beta, binom
from scipy.special import gammaln
from math import exp, log

import csv
# Allow long lines in .cnv files, which can potentially list thousands of SSMs
# in one CNV. According to http://stackoverflow.com/a/15063941, this value can
# be as much as a C long. C longs are guaranteed to accommodate at least this
# value, which is a signed 32-bit int.
csv.field_size_limit(2147483647)

from data import Datum

from tssb import *

def log_factorial(n):
	return gammaln(n + 1)

def log_bin_coeff(n, k):
	return log_factorial(n) - log_factorial(k) - log_factorial(n - k)

def log_binomial_likelihood(x, n, mu):
	return x * log(mu) + (n - x) * log(1 - mu)

def log_beta(a, b):
	return gammaln(a) + gammaln(b) - gammaln(a + b)

def logsumexp(X, axis=None):
    maxes = numpy.max(X, axis=axis)
    return numpy.log(numpy.sum(numpy.exp(X - maxes), axis=axis)) + maxes

def load_data(fname1,fname2):
	# load ssm data
	reader = csv.DictReader(open(fname1,'rU'), delimiter='\t')
	data = dict()  
	for row in reader:
		name = row['gene']
		id = row['id']
		a = [int(x) for x in row['a'].split(',')]
		d = [int(x) for x in row['d'].split(',')]

		mu_r=mu_v=0
		if 'mu_r' in row.keys():
			mu_r = float(row['mu_r'])
			mu_v = float(row['mu_v'])

		data[id] = Datum(name, id, a, d, mu_r, mu_v)
	
	n_ssms = len(data.keys())
	n_cnvs = 0
	
	# load cnv data
	reader = csv.DictReader(open(fname2,'rU'), delimiter='\t')

	for row in reader:
		name=row['cnv']
		id = row['cnv']
		a = [int(x) for x in row['a'].split(',')]
		d = [int(x) for x in row['d'].split(',')]

		data[id] = Datum(name, id, a, d,0.999,0.5)

		ssms = row['ssms']
		if ssms is None: continue
		if len(ssms)>0:
			for ssm in ssms.split(';'):
				tok = ssm.split(',')
				data[tok[0]].cnv.append((data[id],int(tok[1]),int(tok[2])))

	n_cnvs = len(data.keys())-n_ssms
		
	return [data[key] for key in data.keys()], n_ssms, n_cnvs
	
#################################################
## some useful functions to get some info about,
## the tree, used by CNV related computations
def set_node_height(tssb):
	tssb.root['node'].ht=0
	def descend(root,ht):
		for child in root.children():
			child.ht=ht
			descend(child,ht+1)
	descend(tssb.root['node'],1)
	
def set_path_from_root_to_node(tssb):
	nodes = tssb.get_nodes()
	for node in nodes: node.path = node.get_ancestors()

def map_datum_to_node(tssb):
	nodes = tssb.get_nodes()
	for node in nodes:
		for datum in node.get_data():
			datum.node=node
#################################################

def check_bounds(p,l=0.0001,u=.9999):
	if p < l: p=l
	if p > u: p=u
	return p

# removes the empty nodes from the tssb tree
# Does not removes root as it is not required
# root: root of the current tree
# parent: parent of the root
# Note this funciton modifies the sticks so they remain valid.
def remove_empty_nodes(root, parent = None):
	for child in list(root['children']):
		remove_empty_nodes(child, root)
	if (root['node'].get_data() == []):
		if (root['children'] == []): # leaf
			if (parent != None):
				ind = parent['children'].index(root)
				parent['children'].remove(root)
				root['node'].kill()
				parent['sticks'] = delete(parent['sticks'],ind,0)
			return
		else:
			if (parent != None):
				parent_ = root['node'].parent()
				ind = parent['children'].index(root)
				for i,child in enumerate(list(root['children'])):
					parent['children'].append(child)
					toappend = zeros((1,1))
					toappend[0] = root['sticks'][i]
					parent['sticks'] = append(parent['sticks'],toappend,0)
					root['children'].remove(child)
				for child in list(root['node'].children()):
					child._parent = parent_
					parent_.add_child(child)
					root['node'].remove_child(child)
				parent['children'].remove(root)
				parent['sticks'] = delete(parent['sticks'],ind,0)
				root['node'].kill()


def rm_safely(filename):
	try:
	    os.remove(filename)
	except OSError as e:
	    if e.errno == 2: # Ignore "no such file" errors
		pass
	    else:
		raise e

class CorruptZipFileError(Exception):
    pass

class BackupManager(object):
    def __init__(self, filenames):
	self._filenames = filenames
	self._backup_filenames = [os.path.realpath(fn) + '.backup' for fn in self._filenames]

    def save_backup(self):
	for fn, backup_fn in zip(self._filenames, self._backup_filenames):
	    shutil.copy2(fn, backup_fn)

    def restore_backup(self):
	for fn, backup_fn in zip(self._filenames, self._backup_filenames):
	    shutil.copy2(backup_fn, fn)

class StateManager(object):
    default_last_state_fn = 'state.last.pickle'
    default_initial_state_fn = 'state.initial.pickle'

    def __init__(self):
	self._initial_state_fn = StateManager.default_initial_state_fn
	self._last_state_fn = StateManager.default_last_state_fn

    def _write_state(self, state, state_fn):
	with open(state_fn, 'w') as state_file:
	    pickle.dump(state, state_file, protocol=pickle.HIGHEST_PROTOCOL)

    def write_state(self, state):
	self._write_state(state, self._last_state_fn)

    def load_state(self):
	with open(self._last_state_fn) as state_file:
	    return pickle.load(state_file)

    def load_initial_state(self):
	with open(self._initial_state_fn) as state_file:
	    return pickle.load(state_file)

    def write_initial_state(self, state):
	self._write_state(state, self._initial_state_fn)

    def delete_state_file(self):
	rm_safely(self._last_state_fn)

    def state_exists(self):
	return os.path.isfile(self._last_state_fn)


class TreeWriter(object):
    default_archive_fn = 'trees.zip'

    def __init__(self, resume_run = False):
	self._archive_fn = TreeWriter.default_archive_fn
	if resume_run:
	    self._ensure_archive_is_valid()
	else:
	    # Remove file to avoid unwanted behaviour. By the zipfile module's
	    # behaviour, given that we open the file with the "a" flag, if a
	    # non-zip file exists at this path, a zip file will be appended to
	    # the file; otherwise, if the file is already a zip, additional
	    # files will be written into the zip. On a new run, neither case is
	    # something we want.
	    rm_safely(self._archive_fn)

    def _ensure_archive_is_valid(self):
	with zipfile.ZipFile(self._archive_fn) as zipf:
	    if zipf.testzip() is not None:
		raise CorruptZipFileError('Corrupt zip file: %s' % self._archive_fn)

    def _open_archive(self):
	self._archive = zipfile.ZipFile(self._archive_fn, 'a', compression=zipfile.ZIP_DEFLATED, allowZip64=True)

    def _close_archive(self):
	self._archive.close()

    def write_trees(self, trees):
	self._open_archive()
	for tree, idx, llh in trees:
	    is_burnin = idx < 0
	    prefix = is_burnin and 'burnin' or 'tree'
	    treefn = '%s_%s_%s' % (prefix, idx, llh)
	    serialized = pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL)
	    self._archive.writestr(treefn, serialized)
	self._close_archive()

class TreeReader(object):
    def __init__(self, archive_fn):
	self._archive = zipfile.ZipFile(archive_fn)
	infolist = self._archive.infolist()
	tree_info = [t for t in infolist if t.filename.startswith('tree_')]
	burnin_info = [t for t in infolist if t.filename.startswith('burnin_')]

	# Sort by index
	tree_info.sort(key = lambda tinfo: self._extract_metadata(tinfo)[0])
	burnin_info.sort(key = lambda tinfo: self._extract_burnin_idx(tinfo))

	self._trees = []
	self._burnin_trees = []

	for info in tree_info:
	    idx, llh = self._extract_metadata(info)
	    assert idx == len(self._trees)
	    self._trees.append((idx, llh, info))
	for info in burnin_info:
	    idx = self._extract_burnin_idx(info)
	    assert len(burnin_info) + idx == len(self._burnin_trees)
	    self._burnin_trees.append((idx, info))

    def num_trees(self):
	return len(self._trees)

    def close(self):
	self._archive.close()

    def _extract_metadata(self, zinfo):
	tokens = zinfo.filename.split('_')
	idx = int(tokens[1])
	llh = float(tokens[2])
	return (idx, llh)

    def _extract_burnin_idx(self, zinfo):
	idx = int(zinfo.filename.split('_')[1])
	return idx

    def _parse_tree(self, zinfo, remove_empty_vertices=False):
	pickled = self._archive.read(zinfo)
	tree = pickle.loads(pickled)
	if remove_empty_vertices:
	    remove_empty_nodes(tree.root)
	return tree

    def load_tree(self, idx, remove_empty_vertices=False):
	tidx, llh, zinfo = self._trees[idx]
	assert tidx == idx
	return self._parse_tree(zinfo, remove_empty_vertices)

    def load_trees(self, num_trees=None, remove_empty_vertices=False):
	for idx, llh, tree in self.load_trees_and_metadata(num_trees, remove_empty_vertices):
	    yield tree

    def load_trees_and_burnin(self, remove_empty_vertices=False):
	for tidx, zinfo in self._burnin_trees:
	    tree = self._parse_tree(zinfo, remove_empty_vertices)
	    yield (tidx, tree)
	for tidx, llh, zinfo in self._trees:
	    tree = self._parse_tree(zinfo, remove_empty_vertices)
	    yield (tidx, tree)

    def load_trees_and_metadata(self, num_trees=None, remove_empty_vertices=False):
	# Sort by LLH
	trees = sorted(self._trees, key = lambda (tidx, llh, zinfo): llh, reverse=True)

	if num_trees is not None:
	    num_trees = min(num_trees, len(trees))
	    trees = trees[:num_trees]

	for tidx, llh, zinfo in trees:
	    tree = self._parse_tree(zinfo, remove_empty_vertices)
	    yield (tidx, llh, tree)
