from ete2 import *
from util2 import *
from porder import *
import string

count=-1

# writes the latex code to print tree
# fdir: directory where trees are stored
# fout: output file for latex
# k: number of trees to be printed
# for more than 1 tree, filenames are made by appending numbers to fout
def print_trees(fdir,fout,k):
	flist = os.listdir(fdir); 
	if sum([flist[i]=='.DS_Store' for i in arange(len(flist))]):
		flist.remove('.DS_Store')
	k=int(k)
	if (len(flist) < k):
		k=len(flist)
	flist = sort(array(flist,float))[::-1][:k]
	outfile = list(fout)[:-4]
	outfile = string.join(outfile,'')
	i=0
	for fname in flist:
		i=i+1
		fin = './'+fdir+'/'+str(fname)
		fh = open(fin)
		tssb = cPickle.load(fh)
		fh.close()
		print_tree_latex(tssb,outfile+str(i)+'.tex')

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
		tree_file+='-- & '
	else:
		tree_file=tree_file[:-1]
		tree_file+=' & '
	for i in range(len(root.params)):
		tree_file+='{0} & '.format(str(around(root.params[i],3)))
	tree_file=tree_file[:-2]
	tree_file+='\\\\\n'
	for child in root.children():
		tree_file=print_index(child, tree_file)
	return tree_file

# writes the latex code
# tssb: tssb structure of the tree
# fout: output file for latex
def print_tree_latex(tssb,fout):
	global count
	remove_empty_nodes(tssb.root, None)

	fout = open(fout,'w')
	count=-1
	tree_file='\documentclass{article}\n'
	tree_file+='\usepackage{tikz}\n'
	tree_file+='\usepackage{multicol}\n'
	tree_file+='\usetikzlibrary{fit,positioning}\n'
	tree_file+='\\begin{document}\n'
	tree_file+='\\begin{tikzpicture}\n'
	tree_file+='[grow=south, ->, level distance=20mm,\
	every node/.style={circle, minimum size = 8mm, thick, draw =black,inner sep=2mm},\
	every label/.append style={shape=rectangle, yshift=-1mm},\
	level 2/.style={sibling distance=50mm},\
	level 3/.style={sibling distance=20mm},\
	level 4/.style={sibling distance=20mm},\
	every edge/.style={-latex, thick}]\n'
	tree_file+='\n\\'
	tree_file=write_tree(tssb.root['node'], tree_file)
	tree_file+=';\n'
	count=-1
	tree_file+='\end{tikzpicture}\n'
	tree_file+='\\vspace{1cm}\n'	
	tree_file+='\\begin{center}\n'
	tree_file+='\\begin{tabular}{|c|p{5cm}|'
	for i in range(len(tssb.root['node'].params)):
		tree_file+='l|'
	tree_file+='}\n'
	tree_file+='\\hline\n'
	tree_file+='Node & \multicolumn{{1}}{{|c|}}{{SNVs}} & \multicolumn{{{0}}}{{|c|}}{{$\phi$}}\\\\\n'.format(len(tssb.root['node'].params))
	tree_file+='\\hline\n'
	tree_file=print_index(tssb.root['node'], tree_file)
	tree_file+='\\hline\n'
	tree_file+='\\end{tabular}\n'
	tree_file+='\\end{center}\n'
	tree_file+='\end{document}\n'
	fout.write(tree_file)
	fout.close()
	
if __name__ == "__main__":	
	print_trees(sys.argv[1], sys.argv[2], sys.argv[3])
	#python printo_latex.py 'trees' 'latex_output.tex' '1'
