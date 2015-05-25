import numpy
from numpy import *
import cPickle as pickle
import zipfile

import scipy.stats as stat
from scipy.stats import beta, binom
from scipy.special import gammaln
from math import exp, log

import csv

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
	try:
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

	except Exception as e:
		pass
		
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
	wts, nodes = tssb.get_mixture()
	for node in nodes: node.path = node.get_ancestors()

def map_datum_to_node(tssb):
	wts, nodes = tssb.get_mixture()
	for node in nodes:
		for datum in node.get_data():
			datum.node=node
#################################################

	

def check_bounds(p,l=0.0001,u=.9999):
	if p < l: p=l
	if p > u: p=u
	return p

class TreeWriter(object):
    def __init__(self, archive_fn):
	self._archive = zipfile.ZipFile(archive_fn, 'w', compression=zipfile.ZIP_DEFLATED, allowZip64=True)
	self._counter = 0

    def _write_tree(self, tree, tree_fn):
	serialized = pickle.dumps(tree, protocol=pickle.HIGHEST_PROTOCOL)
	self._archive.writestr(tree_fn, serialized)

    def write_tree(self, tree, llh):
	self._write_tree(tree, 'tree_%s_%s' % (self._counter, llh))
	self._counter += 1

    def write_burnin_tree(self, burnin_tree, idx):
	self._write_tree(burnin_tree, 'burnin_%s' % idx)

    def close(self):
	self._archive.close()

class TreeReader(object):
    def __init__(self, archive_fn):
	self._archive = zipfile.ZipFile(archive_fn)
	infolist = self._archive.infolist()
	tree_info = [t for t in infolist if t.filename.startswith('tree_')]

	# Sort by index
	tree_info.sort(key = lambda tinfo: self._extract_metadata(tinfo)[0])
	self._trees = []
	for info in tree_info:
	    idx, llh = self._extract_metadata(info)
	    assert idx == len(self._trees)
	    self._trees.append((idx, llh, info))

    def num_trees(self):
	return len(self._trees)

    def close(self):
	self._archive.close()

    def _extract_metadata(self, zinfo):
	tokens = zinfo.filename.split('_')
	idx = int(tokens[1])
	llh = float(tokens[2])
	return (idx, llh)

    def load_tree(self, idx):
	tidx, llh, zinfo = self._trees[idx]
	assert tidx == idx
	pickled = self._archive.read(zinfo)
	return pickle.loads(pickled)

    def load_trees(self, num_trees=None):
	for idx, llh, tree in self.load_trees_and_llhs(archive_fn, num_trees):
	    yield tree

    def load_trees_and_metadata(self, num_trees=None):
	# Sort by LLH
	trees = sorted(self._trees, key = lambda (tidx, llh, zinfo): llh, reverse=True)

	if num_trees is not None:
	    num_trees = min(num_trees,len(trees))
	    trees = trees[:num_trees]

	for tidx, llh, zinfo in trees:
	    pickled = self._archive.read(zinfo)
	    tree = pickle.loads(pickled)
	    yield (tidx, llh, tree)
