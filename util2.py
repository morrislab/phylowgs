import numpy
from numpy import *

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
