from numpy	  import *
from numpy.random import *
from tssb		import *
from util		import *
from util2 import *

from ete2 import *

from subprocess import call
import sys

ctr=0
def print_top_trees(tree_archive,fout,k=5):
	global ctr;
	fout = open(fout,'w')
	tree_reader = TreeReader(tree_archive)

	for idx, (tidx, llh, tree) in enumerate(tree_reader.load_trees_and_metadata(k)):
			ctr=0
			remove_empty_nodes(tree.root, None)
			# print top K trees in ascii
			print_best_tree(tree,fout)
			
	tree_reader.close()
	fout.close()	

def print_best_tree(tssb,fout):
	nodes = tssb.get_nodes()
	nnodes = sum( [ 1 for node in nodes if len(node.get_data()) ] )
	
	#fout=open(fout,'w')
	t = Tree();t.name='0'
	fout.write('id, \t phi, \t nChildren, \t nGenes, \t genes \n')
	print_node2(tssb.root,None,t,fout)
	fout.write('\n\n')
	fout.write(t.get_ascii(show_internal=True))
	fout.write('\n\n')	
	fout.write('Number of non-empty nodes in the tree: ' +repr(nnodes))	
	fout.write('\n\n\n')	
	#fout.close()

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
	out_str = str(node_name) + ',\t' + str(around(node['node'].params,3)) + ',\t' + str(len(node['children'])) + ',\t' + str(len(genes)) + ',\t' + gnames +  '\n'
	fout.write(out_str)
	
	for child in node['children']:
		name_string = str(ctr)#+'('+str(len(child['node'].get_data()))+')'
		print_node2(child, node_name,tree.add_child(name=name_string),fout)
	
### printing stuff #################
def print_best_tree_pdf(tssb,fout,score=0):
	#wts, nodes = tssb.get_mixture()
	#w = dict([(n[1], n[0]) for n in zip(wts,nodes)])
	print_tree_latex(tssb,fout,score)	


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
	tree_file+=' & '
	for i in range(len(root.params)):
		tree_file+='{0}, '.format(str(around(root.params[i],3)))
	tree_file=tree_file[:-2]
	tree_file+='\\\\\n'
	for child in root.children():
		tree_file=print_index(child, tree_file)
	return tree_file

# writes the latex code
# tssb: tssb structure of the tree
# fout: output file for latex
def print_tree_latex(tssb,fout,score):
	global count

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
	tree_file+='\\begin{tabular}{|c|p{5cm}|p{5cm}|'
	for i in range(len(tssb.root['node'].params)):
		tree_file+='l|'
	tree_file+='}\n'
	tree_file+='\\hline\n'
	tree_file+='Node & \multicolumn{{1}}{{|c|}}{{Mutations}} & \multicolumn{{1}}{{|c|}}{{Frequencies}} \\\\\n'.format(len(tssb.root['node'].params))
	tree_file+='\\hline\n'
	tree_file=print_index(tssb.root['node'], tree_file)
	tree_file+='\\hline\n'
	tree_file+='\\end{tabular}\n'
	tree_file+='};\n'
	tree_file+='\\end{tikzpicture}\n'
	tree_file+='};\n'	
	#tree_file+='\\node at (b.south) [anchor=north,yshift=-.5cm]{Posterior probability: ' + str(score) + '};\n'
	tree_file+='\\end{tikzpicture}\n'
	tree_file+='\end{document}\n'
	fout.write(tree_file)
	fout.close()	

