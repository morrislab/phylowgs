import cPickle

from numpy	  import *
from numpy.random import *
from tssb		import *
from util2 import *

from ete2 import *
import heapq

from subprocess import call

def compute_lineages(fdir,fin1,fin2,fout):
	codes, n_ssms, n_cnvs = load_data(fin1,fin2)	
	m = len(codes) # number of SSMs+CNVs
	
	flist = os.listdir(fdir)
	if sum([flist[i]=='.DS_Store' for i in arange(len(flist))]): flist.remove('.DS_Store')
	ns = len(flist) #number of MCMC samples
	
	G = dict()
	datalist = dict([(datum.name,str(i)) for i,datum in enumerate(codes)])
	
	for idx,fname in enumerate(flist):
		#if idx>100:continue
		f=open('./'+fdir+'/'+str(fname))
		tssb = cPickle.load(f)
		f.close()
	
		wts,nodes = tssb.get_mixture()
		genotype = dict([(id(node),'') for i,node in enumerate(nodes)])
		gtypelist = []
		# process descendants
		def descend(root):			
			data = root.get_data()
			ndata = len(data)
			gtype=''
			if ndata>0:
				if root.parent() is not None:
					gtype = genotype[id(root.parent())]#??
				for datum in data: 
						gtype += datalist[datum.name]+'_'
				genotype[id(root)] = gtype
				gtypelist.append(_sort(gtype))
			for child in root.children():
				descend(child)
		descend(tssb.root['node'])

		sgtypelist = sort_and_merge(gtypelist)
		if sgtypelist in G: 
			G[sgtypelist].append(idx)
		else:
			G[sgtypelist] = [idx]


	post_trees = []
	flist = os.listdir(fdir)
	if sum([flist[i]=='.DS_Store' for i in arange(len(flist))]): flist.remove('.DS_Store')
	for ptree in G.keys():
		tssb_list = G[ptree]
		prob =round(len(tssb_list)*1./ns,4)
		#print 'posterior probability: ' + repr(prob)
		idx = tssb_list[0]
		heapq.heappush(post_trees,(-1.*prob,idx))


	# print the trees in ascii format
	'''
	fout = open(fout,'w')
	while(len(post_trees)):
		score,idx = heapq.heappop(post_trees)
		print_best_tree('./'+fdir+'/'+str(flist[idx]),fout,-1.*score)
	fout.flush()
	fout.close()
	'''

	# print the trees in latext format
	fout = open(fout,'w')
	fidx=0
	while(len(post_trees)):
		score,idx = heapq.heappop(post_trees)
		print_best_tree('./'+fdir+'/'+str(flist[idx]),'./latex/'+str(fidx)+'.tex',-1.*score)
		

		# system call pdflatex?
		call(['pdflatex', './latex/'+str(fidx)+'.tex','-output-directory=./latex/'])

		fidx+=1

	
	fout.flush()
	fout.close()




def _sort(str):
	str = sort(str.split('_'))
	sstr=''
	for s in str: sstr+= s+'_'
	return sstr	
def sort_and_merge(gtypelist):
	sglist = sort(gtypelist)
	sstr=''
	for s in sglist: sstr+= s+';'
	return sstr	


### printing stuff #################

def print_best_tree(fin,fout,score):
	fh = open(fin)
	tssb = cPickle.load(fh)
	fh.close()
	
	remove_empty_nodes(tssb.root, None) # removes empty leaves
	
	#wts, nodes = tssb.get_mixture()
	#w = dict([(n[1], n[0]) for n in zip(wts,nodes)])
	
	print_tree_latex(tssb,fout,score)

ctr=0
def print_best_tree1(fin,fout,score):
	global ctr
	ctr=0
	fh = open(fin)
	tssb = cPickle.load(fh)
	fh.close()
	
	remove_empty_nodes(tssb.root, None) # removes empty leaves
	
	#wts, nodes = tssb.get_mixture()
	#w = dict([(n[1], n[0]) for n in zip(wts,nodes)])
	
	t = Tree();t.name='0'
	fout.write('id, \t SSMs \n')
	print_node2(tssb.root,None,t,fout)
	fout.write('\n\n')
	fout.write(t.get_ascii(show_internal=True))
	fout.write('\n\n')	
	fout.write('Posterior probability: ' + repr(score))
	fout.write('\n\n')	
	fout.write('###############################################################\n\n')		


	

def print_node2(node, parent,tree,fout):
	global ctr;
	num_data = node['node'].num_data()
	node_name  = ctr ; ctr+=1;
	
	genes = node['node'].get_data()
	gnames = ''
	if len(genes)>0:
		gnames = genes[0].id#name
		for g in arange(1,len(genes)):
			gnames = gnames + '; ' + genes[g].id#name;
	out_str = str(node_name) + ',\t' +  gnames +  '\n'
	fout.write(out_str)
	
	for child in node['children']:
		name_string = str(ctr)#+'('+str(len(child['node'].get_data()))+')'
		print_node2(child, node_name,tree.add_child(name=name_string),fout)
		

# removes the empty nodes from the tssb tree
# Does not removes root as it is not required
# root: root of the current tree
# parent: parent of the root
def remove_empty_nodes(root, parent):
	for child in list(root['children']):
		remove_empty_nodes(child, root)
	if (root['node'].get_data() == []):
		if (root['children'] == []): # leaf
			if (parent != None):
				parent['children'].remove(root)
				root['node'].kill()
			return
		else:
			if (parent != None):
				parent_ = root['node'].parent()
				for child in list(root['children']):
					parent['children'].append(child)
					root['children'].remove(child)
				for child in list(root['node'].children()):
					child._parent = parent_
					parent_.add_child(child)
					root['node'].remove_child(child)
				parent['children'].remove(root)
				root['node'].kill()



################ LATEX PRINTING ######################
global count
# writes code for tree
# root: root of the tree
# tree_file: string with latex code
def write_tree(root, tree_file):
	global count
	count+=1
	tree_file+='node {{{0}}}'.format(count)
	for child in root.children():
		tree_file+='child {'
		tree_file=write_tree(child, tree_file)
		tree_file+='}'
	return tree_file

# writes code for index
# root: root of the tree
# tree_file: string with latex code
def print_index(root, tree_file):
	global count
	count+=1
	tree_file+='{0} & '.format(count)
	ssms=''
	for datum in root.get_data():
		ssms+='{0}, '.format(datum.name)
	tree_file+=ssms.strip().strip(',')
	if root.get_data()==[]:
		tree_file+='-- '
	#else:
	#	tree_file=tree_file[:-1]
		#tree_file+=' & '
	#for i in range(len(root.params)):
	#	tree_file+='{0} & '.format(str(around(root.params[i],3)))
	#tree_file=tree_file[:-2]
	tree_file+='\\\\\n'
	for child in root.children():
		tree_file=print_index(child, tree_file)
	return tree_file

# writes the latex code
# tssb: tssb structure of the tree
# fout: output file for latex
def print_tree_latex(tssb,fout,score):
	global count
	#remove_empty_nodes(tssb.root, None)

	fout = open(fout,'w')
	count=-1
	#tree_file='\documentclass{article}\n'
	tree_file='\documentclass{standalone}\n'	
	tree_file+='\usepackage{tikz}\n'
	tree_file+='\usepackage{multicol}\n'
	tree_file+='\usetikzlibrary{fit,positioning}\n'
	tree_file+='\\begin{document}\n'
	tree_file+='\\begin{tikzpicture}\n'
	tree_file+='\\node (a) at (0,0){\n'
	tree_file+='\\begin{tikzpicture}\n'	
	tree_file+='[grow=east, ->, level distance=20mm,\
	every node/.style={circle, minimum size = 8mm, thick, draw =black,inner sep=2mm},\
	every label/.append style={shape=rectangle, yshift=-1mm},\
	level 2/.style={sibling distance=50mm},\
	level 3/.style={sibling distance=20mm},\
	level 4/.style={sibling distance=20mm},\
	every edge/.style={-latex, thick}]\n'
	tree_file+='\n\\'
	tree_file=write_tree(tssb.root['node'], tree_file)
	tree_file+=';\n'
	tree_file+='\\end{tikzpicture}\n'
	tree_file+='};\n'	
	count=-1
	tree_file+='\\node (b) at (a.south)[anchor=north,yshift=-.5cm]{\n'
	tree_file+='\\begin{tikzpicture}\n'	
	tree_file+='\\node (table){\n'
	tree_file+='\\begin{tabular}{|c|p{5cm}|'
	for i in range(len(tssb.root['node'].params)):
		tree_file+='l|'
	tree_file+='}\n'
	tree_file+='\\hline\n'
	tree_file+='Node & \multicolumn{{1}}{{|c|}}{{Mutations}}\\\\\n'.format(len(tssb.root['node'].params))
	tree_file+='\\hline\n'
	tree_file=print_index(tssb.root['node'], tree_file)
	tree_file+='\\hline\n'
	tree_file+='\\end{tabular}\n'
	tree_file+='};\n'
	tree_file+='\\end{tikzpicture}\n'
	tree_file+='};\n'	
	tree_file+='\\node at (b.south) [anchor=north,yshift=-.5cm]{Posterior probability: ' + str(score) + '};\n'
	tree_file+='\\end{tikzpicture}\n'
	tree_file+='\end{document}\n'
	fout.write(tree_file)
	fout.close()	
		
if __name__ == "__main__":


	dir = './best/' # the folder name with best trees	

	fin1='ssm_data.txt'
	fin2='cnv_data.txt'
	
	compute_lineages(dir,fin1,fin2,'postk')	
