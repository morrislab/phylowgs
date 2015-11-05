###### code to sample from the paramater posterior p(\phi | data) ########

import numpy
from numpy import *
from data import Datum

from tssb import *

from util import dirichletpdfln
from numpy.random import dirichlet

import subprocess as sp

import util2 as u2
import os

def get_c_fnames(tmp_dir):
	def _make_c_fname(name):
		fname = 'c_%s.txt' % (name)
		return os.path.join(tmp_dir, fname)

	FNAME_C_TREE = _make_c_fname('tree')
	FNAME_C_DATA_STATES = _make_c_fname('data_states')
	FNAME_C_PARAMS = _make_c_fname('params')
	FNAME_C_MH_ARATIO = _make_c_fname('mh_ar')

	return (FNAME_C_TREE, FNAME_C_DATA_STATES, FNAME_C_PARAMS, FNAME_C_MH_ARATIO)

# done for multi-sample
def metropolis(tssb,iters=1000,std=0.01,burnin=0,n_ssms=0,n_cnvs=0,fin1='',fin2='',rseed=1, ntps=5, tmp_dir='.'):
	wts, nodes = tssb.get_mixture()

	# file names
	FNAME_SSM_DATA = fin1
	FNAME_CNV_DATA = fin2
	NTPS = str(ntps)
	FNAME_C_TREE, FNAME_C_DATA_STATES, FNAME_C_PARAMS, FNAME_C_MH_ARATIO = get_c_fnames(tmp_dir)

	## initialize the MH sampler###########
	#for tp in arange(ntps): 
	#	sample_cons_params(tssb,tp)
	#	update_params(tssb,tp)
	######################################
	
	## prepare to call the c++ code ###########
	u2.set_node_height(tssb)
	write_tree(tssb,n_ssms,FNAME_C_TREE) #write the current tree to the disk
	u2.map_datum_to_node(tssb)
	
	write_data_state(tssb,FNAME_C_DATA_STATES) # this is need for binomial parameter computations
	###########################################
	
	MH_ITR = str(iters)
	MH_STD = str(std)
	N_SSM_DATA = str(n_ssms)
	N_CNV_DATA = str(n_cnvs)
	NNODES = str(len(nodes))
	TREE_HEIGHT = str(max([node.ht for node in nodes])+1)
	
	script_dir = os.path.dirname(os.path.realpath(__file__))
	sp.check_call(['%s/mh.o' % script_dir, MH_ITR, MH_STD, N_SSM_DATA, N_CNV_DATA, NNODES, TREE_HEIGHT, FNAME_SSM_DATA, FNAME_CNV_DATA, FNAME_C_TREE, FNAME_C_DATA_STATES, FNAME_C_PARAMS,FNAME_C_MH_ARATIO, NTPS])
	ar = str(loadtxt(FNAME_C_MH_ARATIO,dtype='string'))
	update_tree_params(tssb,FNAME_C_PARAMS) # update the tree with the new parameters sampled using the c++ code
	
	return ar

# done for multi-sample
def write_tree(tssb,n_ssms,fname):
	fh=open(fname,'w')
	wts,nodes=tssb.get_mixture()
	did_int_dict=dict()
	for dat in tssb.data:
		if dat.id[0]=='s':
			did_int_dict[dat.id]=int(dat.id[1:])
		else:
			did_int_dict[dat.id]=n_ssms+int(dat.id[1:])
	
	def descend(root):		
		for child in root.children():			
			descend(child)
		
		# write data#
		cids=''
		for child in root.children():cids+=str(child.id)+','
		cids=cids.strip(',')
		if cids=='': cids=str(-1)
		
		dids=''
		for dat in root.get_data():dids+=str(did_int_dict[dat.id])+','
		dids=dids.strip(',')
		if dids=='': dids=str(-1)
		
		line = str(root.id) + '\t' + list_to_string(root.params) + '\t' + list_to_string(root.pi) + '\t' + str(len(root.children())) + '\t'  + cids + '\t' + str(len(root.get_data())) + '\t' + dids + '\t' +  str(root.ht)
		fh.write(line)
		fh.write('\n')
		fh.flush()
		###############
	
	descend(tssb.root['node'])
	fh.flush()

	fh.close()


def list_to_string(p):
	o=''
	for pp in p:o+=str(pp)+','
	return o.strip(',')
	
# no changes for multi-sample
# data/node state format (parameter independent dot-product weights)
# datum_id	node_id_1,pi,nr,nv;node_id_2,pi,nr,nv;....	
# these weights are used to compute data log-likelihood			
def write_data_state(tssb,fname):
	fh = open(fname,'w')
	wts,nodes=tssb.get_mixture()
	
	for dat in tssb.data:
		if not dat.cnv: continue # nothing to do for CNVs
		if not dat.node: continue # todo: this won't happen
		poss_n_genomes = dat.compute_n_genomes(0)
		if poss_n_genomes[0][1] == 0:
			nv = (False,True)
		elif poss_n_genomes[1][1] == 0:
			nv = (True,False)
		else:
			nv = (True,True)
		for node in nodes:
		
			ssm_node = node.path[-1]
			mr_cnv = find_most_recent_cnv(dat,node)
			ancestors = node.get_ancestors()
            
			dat.state1 = '' # maternal
			dat.state2 = '' # paternal
			if (not ssm_node in ancestors) and (not mr_cnv):
				dat.state1 += str(node.id) + ',' + str(2) + ',' + str(0) + ';'
				dat.state2=dat.state1
			elif ssm_node in ancestors and (not mr_cnv):
				dat.state1 += str(node.id) + ',' + str(1) + ',' + str(1) + ';'
				dat.state2=dat.state1
			elif (not ssm_node in ancestors) and mr_cnv:
				dat.state1 += str(node.id) + ',' + str(mr_cnv[1]+mr_cnv[2]) + ',' + str(0) + ';'
				dat.state2=dat.state1
			elif ssm_node in ancestors and mr_cnv:
				if ssm_node in mr_cnv[0].node.get_ancestors():
					if nv == (False,True):
						dat.state2 += str(node.id) + ',' + str(mr_cnv[2]) + ',' + str(mr_cnv[1]) + ';' # paternal
						dat.state1=dat.state2
					elif nv == (True, False):
						dat.state1 += str(node.id) + ',' + str(mr_cnv[1]) + ',' + str(mr_cnv[2]) + ';' # maternal
						dat.state2 = dat.state1
					else:
						dat.state1 += str(node.id) + ',' + str(mr_cnv[1]) + ',' + str(mr_cnv[2]) + ';' # maternal
						dat.state2 += str(node.id) + ',' + str(mr_cnv[2]) + ',' + str(mr_cnv[1]) + ';' # paternal
					
				else:
					dat.state1 += str(node.id) + ',' + str(max(0,mr_cnv[1]+mr_cnv[2]-1)) + ',' + str(min(1,mr_cnv[1]+mr_cnv[2])) + ';'
					dat.state2 = dat.state1 
			else:
				print "PANIC"
			
			fh.write(str(dat.id[1:]) + '\t' + dat.state1.strip(';') + '\t' + dat.state2.strip(';'))
			fh.write('\n')
		
	fh.flush()
	fh.close()

# done for multi-sample	
def find_most_recent_cnv(dat,nd):
	out = None
	for n in nd.get_ancestors()[::-1]:
		if n in [x[0].node for x in dat.cnv]:
			out = [x for x in dat.cnv if x[0].node == n][0]
			break
	return out

# done for multi sample
def update_tree_params(tssb,fname):
	wts, nodes = tssb.get_mixture()
	ndict = dict()
	for node in nodes: ndict[node.id]=node
	
	fh=open(fname)
	params=[line.split() for line in fh.readlines()]
	fh.close()
	
	for p in params:
		ndict[int(p[0])].params = string_to_list(p[1])
		ndict[int(p[0])].pi = string_to_list(p[2])
	#params=loadtxt('c_params.txt')
	#for p in params:
	#	ndict[p[0]].params = p[1]
	#	ndict[p[0]].pi = p[2]
	

def string_to_list(p):
	p=p.strip(',')
	return array([float(pp) for pp in p.split(',')])

# done for multi-sample
# tree-structured finite-dimensional stick breaking
def sample_cons_params(tssb,tp):
	def descend(root,tp):

		if root.parent() is None:
			root.params1[tp] = 1
			root.pi1[tp] = root.params1[tp]*rand(1) # break nu stick
		r = root.params1[tp]-root.pi1[tp] #mass assigned to children
		p = rand(len(root.children()));p=r*p*1./sum(p)
		index=0
		for child in root.children():			
			child.params1[tp] = p[index]# break psi sticks			
			child.pi1[tp] = child.params1[tp]*(rand(1)**(len(child.children())>0)) # break nu stick
			index+=1
		for child in root.children():			
			descend(child,tp)	

	descend(tssb.root['node'],tp)

# done for multi-sample
def update_params(tssb,tp):
	def descend(root,tp):			
		for child in root.children():
			descend(child,tp)	
		root.params[tp] = root.params1[tp]
		root.pi[tp] = root.pi1[tp]
	descend(tssb.root['node'],tp)
	
	
	
	
	
	
###### old code, not in use #############
# data/node state format (parameter independent dot-product weights)
# datum_id	node_id_1,pi,nr,nv;node_id_2,pi,nr,nv;....	
# these weights are used to compute data log-likelihood			
def write_data_state1111(tssb):
	fh = open('c_data_states.txt','w')
	wts,nodes=tssb.get_mixture()
	
	for dat in tssb.data:
		if not dat.cnv: continue # nothing to do for CNVs
		
		if not dat.node: continue # todo: this won't happen
		
		ancestors = dat.node.get_ancestors() # path from root to ssm node
	
		mr_cnv = dat.cnv[0] # CNV corresponding to the SSM
		
		dat.state1 = '' # maternal
		dat.state2 = '' # paternal
		
		# do this until we encounter the SSM node,
		# i.e., along the path from root to the SSM node
		visited_cnv = False
		for node in ancestors:
		
			if node != mr_cnv[0].node and visited_cnv==False: # until CNV is encountered
				dat.state1 += str(node.id) + ',' + str(2) + ',' + str(0) + ';'
			else:
				visited_cnv = True
				dat.state1 += str(node.id) + ',' + str(mr_cnv[1]+mr_cnv[2]) + ',' + str(0) + ';'
			dat.state2=dat.state1	
		
		# do this after the SSM node, i.e, for all nodes in the subtree below the SSM node
		# [node_id, nr, nv] format
		def descend(n,d):
			if n == mr_cnv[0].node:
				d.state1 += str(n.id) + ',' + str(mr_cnv[1]) + ',' + str(mr_cnv[2]) + ';' # maternal
				d.state2 += str(n.id) + ',' + str(mr_cnv[2]) + ',' + str(mr_cnv[1]) + ';' # paternal
			else:
				d.state1 += str(n.id) + ',' + str(mr_cnv[1]+mr_cnv[2]-1) + ',' + str(1) + ';'
				d.state2 = d.state1
			for child in n.children():
				descend(child,d)
		
		# traverse the tree below the ssm node
		for child in node.children(): descend(child,dat)
		
		fh.write(str(dat.id[1:]) + '\t' + dat.state1.strip(';') + '\t' + dat.state2.strip(';'))
		fh.write('\n')
		
	fh.flush()
	fh.close()





