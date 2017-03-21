#!/usr/bin/env python2
import cPickle

from numpy	  import *
from numpy.random import *
from tssb		import *
from util2 import *

from ete2 import *
import heapq

from subprocess import call

import argparse
import sys

def compute_lineages(archive_fn, num_trees, fin1, fin2):
	codes, n_ssms, n_cnvs = load_data(fin1,fin2)	
	m = len(codes) # number of SSMs+CNVs
	tree_reader = TreeReader(archive_fn)
	ns = tree_reader.num_trees() #number of MCMC samples
	
	G = dict()
	datalist = dict([(datum.name,str(i)) for i,datum in enumerate(codes)])
	
	for idx, llh, tssb in tree_reader.load_trees_and_metadata():
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
	for ptree in G.keys():
		tssb_list = G[ptree]
		prob =round(len(tssb_list)*1./ns,4)
		#print 'posterior probability: ' + repr(prob)
		idx = tssb_list[0]
		# Note that trees aren't ordered by likelihood -- only posterior
		# probabaility.
		heapq.heappush(post_trees,(-1.*prob,tssb_list))


	# print the trees in latex format
	try:
		os.mkdir('posterior_trees')
	except OSError, e:
		if e.errno == 17: # Directory exists
			pass
		else:
			raise e
	fidx=0

	if num_trees is None:
		num_trees = len(post_trees)

	while len(post_trees) > 0 and fidx < num_trees:
		score,idx = heapq.heappop(post_trees)
		score = -score

		# aggregate frequencies
		freqs = dict()
		for id1 in idx:
			tssb = tree_reader.load_tree(id1)
			remove_empty_nodes(tssb.root, None)

			def descend(root):
				for child in root.children():
					descend(child)
				names=''
				for dat in root.get_data():names+=dat.name+';'
				names=names.strip(';')
				if names in freqs:
					freqs[names].append(root.params)
				else:
					if len(names)>0 or root.parent() is None:
						freqs[names]=[]
						freqs[names].append(root.params)
			descend(tssb.root['node'])

		tex_fn = 'posterior_trees/tree_%s_%s.tex' % (fidx, score)
		print_best_tree(tree_reader.load_tree(idx[0]), tex_fn, score, freqs)

		# Call pdflatex. To permit it to find standalone.* files,
		# change into PhyloWGS directory to run the command, then
		# change back to previous directory.
		script_dir = os.path.dirname(os.path.realpath(__file__))
		old_wd = os.getcwd()
		os.chdir(script_dir)
		try:
			call(['pdflatex', '-interaction=nonstopmode', '-output-directory=%s/posterior_trees/' % old_wd, '%s/%s' % (old_wd, tex_fn)])
		except OSError:  # pdflatex not available, do not die
			print >> sys.stderr, 'pdflatex not available'
		os.chdir(old_wd)

		fidx+=1
	tree_reader.close()

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

def print_best_tree(tssb, fout, score, freqs):
	remove_empty_nodes(tssb.root, None) # removes empty leaves
	
	#wts, nodes = tssb.get_mixture()
	#w = dict([(n[1], n[0]) for n in zip(wts,nodes)])
	
	print_tree_latex(tssb,fout,score,freqs)

ctr=0
def print_best_tree1(tssb,fout,score):
	global ctr
	ctr=0
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
def print_index(root, tree_file,freqs):
	global count
	count+=1
	tree_file+='{0} & '.format(count)

	num_ssms = 0
	num_cnvs = 0
	mutations = root.get_data()
	for mut in mutations:
		if mut.id.startswith('s'):
			num_ssms += 1
		elif mut.id.startswith('c'):
			num_cnvs += 1
		else:
			raise Exception('Unknown mutation ID type: %s' % mut.id)

	tree_file += '%s & %s &' % (num_ssms, num_cnvs)

	# print params
	names=''
	for dat in root.get_data():names+=dat.name+';'
	names=names.strip(';')
	freq = array(freqs[names])
	for i in range(len(root.params)):
		#tree_file+='{0} & '.format(str(around(root.params[i],3)))
		tree_file+='{0} & '.format( str(around(mean(freq[:,i]),3)) + ' $\pm$ ' + str(around(std(freq[:,i]),3)) )
	tree_file=tree_file[:-2]
	tree_file+='\\\\\n'

	for child in root.children():
		tree_file=print_index(child, tree_file,freqs)
	return tree_file

# writes the latex code
# tssb: tssb structure of the tree
# fout: output file for latex
def print_tree_latex(tssb,fout,score,freqs):
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
	tree_file+='\\begin{tabular}{|c|l|l|'
	for i in range(len(tssb.root['node'].params)):
		tree_file+='l|'
	tree_file+='}\n'
	tree_file+='\\hline\n'
	tree_file+='Node & \multicolumn{{1}}{{|c|}}{{SSMs}} & \multicolumn{{1}}{{|c|}}{{CNVs}} & \multicolumn{{{0}}}{{|c|}}{{Clonal frequencies}}\\\\\n'.format(len(tssb.root['node'].params))
	tree_file+='\\hline\n'
	tree_file=print_index(tssb.root['node'], tree_file,freqs)
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
	parser = argparse.ArgumentParser(
		description='Plot posterior trees resulting from PhyloWGS run',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument('--num-trees', '-n', dest='num_trees', type=int,
		help='Only output given number of trees')
	parser.add_argument('ssm_file',
		help='File listing SSMs (simple somatic mutations, i.e., single nucleotide variants. For proper format, see README.md.')
	parser.add_argument('cnv_file',
		help='File listing CNVs (copy number variations). For proper format, see README.md.')
	args = parser.parse_args()

	compute_lineages(TreeWriter.default_archive_fn, args.num_trees, args.ssm_file, args.cnv_file)
