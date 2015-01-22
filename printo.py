import cPickle

from numpy	  import *
from numpy.random import *
from tssb		import *
from util		import *
from util2 import *

from ete2 import *

ctr=0
def print_top_trees(fdir,fout,k=5):
	global ctr;
	flist = os.listdir(fdir); 
	if sum([flist[i]=='.DS_Store' for i in arange(len(flist))]): flist.remove('.DS_Store')
	if len(flist)<k:k=len(flist)
	flist = sort(array(flist,float))[::-1][:k]
	fout = open(fout,'w')
	for fname in flist:
			ctr=0
			print_best_tree('./'+fdir+'/'+str(fname),fout)
	fout.close()	

def print_best_tree(fin,fout):
	fh = open(fin)
	tssb = cPickle.load(fh)
	fh.close()
	
	wts, nodes = tssb.get_mixture()
	w = dict([(n[1], n[0]) for n in zip(wts,nodes)])
	nnodes = sum( [ 1 for node in nodes if len(node.get_data()) ] )
	
	#fout=open(fout,'w')
	t = Tree();t.name='0'
	fout.write('id, \t phi, \t nChildren, \t nGenes, \t genes \n')
	print_node2(tssb.root,None,t,w,fout)
	fout.write('\n\n')
	fout.write(t.get_ascii(show_internal=True))
	fout.write('\n\n')	
	fout.write('Number of non-empty nodes in the tree: ' +repr(nnodes))	
	fout.write('\n\n\n')	
	#fout.close()

def print_node2(node, parent,tree,wts,fout):
	global ctr;
	num_data = node['node'].num_data()
	node_name  = ctr ; ctr+=1;
	
	genes = node['node'].get_data()
	gnames = ''
	if len(genes)>0:
		gnames = genes[0].id#name
		for g in arange(1,len(genes)):
			gnames = gnames + '; ' + genes[g].id#name;
	out_str = str(node_name) + ',\t' + str(around(node['node'].params,3)) + ',\t' + str(len(node['children'])) + ',\t' + str(len(genes)) + ',\t' + gnames +  '\n'
	fout.write(out_str)
	
	for child in node['children']:
		name_string = str(ctr)#+'('+str(len(child['node'].get_data()))+')'
		print_node2(child, node_name,tree.add_child(name=name_string),wts,fout)