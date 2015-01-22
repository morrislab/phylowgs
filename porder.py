# code to generate partial order plots

import cPickle

from numpy	  import *
from numpy.random import *
from tssb		import *
from util		import *
from util2 import *

from ete2 import *
import pygraphviz as pgv

from cc import *


# print partial order
# fdir: directory where trees are stored
# dname: input data file name
# fout: output file for graph figure
def print_porder(fdir,dname,fout):
	codes = load_data(dname)
	m = len(codes) # number of genes / instances
	glist = dict([(datum.name,i) for i,datum in enumerate(codes)])
	nparams = dict([(datum.name,0) for i,datum in enumerate(codes)]) # clonal frequencies
	 
	W = zeros((m,m)) # parent-child matrix
	C = zeros((m,m)) #cluster matrix
	S = zeros((m,m)) #similarity matrix (similar to pyclone)
	flist = os.listdir(fdir)
	#if len(where(flist=='.DS_Store')[0]): flist.remove('.DS_Store')
	if sum([flist[i]=='.DS_Store' for i in arange(len(flist))]): flist.remove('.DS_Store')
	flist = sort(array(flist,float))[::-1]
	ns = len(flist) #number of MCMC samples
	
	for idx,fname in enumerate(flist):
		f=open('./'+fdir+'/'+str(fname))
		tssb = cPickle.load(f)
		f.close()
		
		def descend(root):	
			for child in root.children():
				cids = [glist[datum.name] for datum in child.get_data()]
				ancestors = [child.parent()] 
				for anc in ancestors:				
					pids = [glist[datum.name] for datum in anc.get_data()]
					ids=[(p,c) for p in pids for c in cids]
					for id in ids: W[id[0],id[1]]+=1
				descend(child)
		descend(tssb.root['node'])
	
		
		#cluster info
		wts, nodes = tssb.get_mixture()
		for node in nodes:
			data = node.get_data()
			for datum in data: nparams[datum.name]+=around(node.params,2)
			
			pids = [glist[datum.name] for datum in data];npids=len(pids)			
			nids = setdiff1d(arange(m),pids);nnids=len(nids)
			ids=[(pids[i],pids[j]) for i in arange(npids) for j in arange(npids)]
			for id in ids: 
				C[id[0],id[1]]+=1 # same cluster
				S[id[0],id[1]]+=1 # same cluster
			ids=[(i,j) for i in pids for j in nids]
			for id in ids: C[id[0],id[1]]-=1; C[id[1],id[0]]-=1 # different cluster
	
	for key in nparams.keys():nparams[key]=nparams[key]*1./ns
	
	# save similarity matrix (similar to pyclone)
	#S[arange(m),arange(m)]=ns
	#savetxt('sim_matrix',S,fmt='%s')
	
	#correlation clustering
	#C = around(C*1./ns)
	C[arange(m),arange(m)]=ns
	C = around(array(cc_lp(C)))
	z = get_cluster_assignments(C)
	
	W = around(W*1./ns,2)

	#glist = dict([(i,datum.name+' ('+str(round(nparams[datum.name],2)) + ')') for i,datum in enumerate(codes)])
	glist = dict([(i,datum.name) for i,datum in enumerate(codes)])
	
	draw_graph(W,glist,z,fout)	


def get_cluster_assignments(C):
	z = -1*ones(C.shape[0])	
	ctr=0
	while 1:
		id=where(z==-1)[0]
		if len(id)==0: break
		id=id[0]
		z[where(C[id]==1)[0]]=ctr
		ctr+=1
	return array(z,int)
	
def draw_graph(W,glist,z,fout='x'):
	m = W.shape[0]
	G = pgv.AGraph(directed=True)
		
	nodeids = argsort(z)
	
	colors = ['red','green','purple','indigo','cyan','orange','black','brown','blue','purple','maroon']	
	#colors = ['red','green','blue','black','brown','red','green','blue','black','brown','red','green','blue','black','brown']
	
	nodes = glist.keys()
	for i,n in enumerate(nodes):
		G.add_node(glist[n],color=colors[z[n]],penwidth=2)
		
	edges = where(W!=0)
	edges = zip(edges[0],edges[1])
	for e in edges:
		if W[e[0],e[1]] < 0.09: continue
		if W[e[0],e[1]] >= 0: 
			G.add_edge(glist[e[0]],glist[e[1]],penwidth=5*(W[e[0],e[1]]))
		else:
			G.add_edge(glist[e[0]],glist[e[1]],color='gray',penwidth=5*(W[e[0],e[1]]))
	G.layout()
	G.draw(fout,prog='dot') 	


if __name__ == "__main__":	
	print_porder(sys.argv[1], sys.argv[2], sys.argv[3])

# python2.6 porder.py './best' 'data.txt' 'partial_order.pdf'

	
	
	
