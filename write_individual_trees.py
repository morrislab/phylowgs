#!/usr/bin/env python2
import cPickle

from numpy    import *
from numpy.random import *
from tssb    import *
from util2 import *

from ete2 import *
import heapq

from subprocess import call

import argparse

def compute_lineages(fdir,fin1,fin2):
  flist = os.listdir(fdir)
  if sum([flist[i]=='.DS_Store' for i in arange(len(flist))]): flist.remove('.DS_Store')
  
  for idx, fname in enumerate(flist):
    #if idx>100:continue
    with open(fdir+'/'+str(fname)) as f:
      tssb = cPickle.load(f)
      remove_empty_nodes(tssb.root, None)
      print_tree_latex(tssb, 'latex/%s.tex' % fname)

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

  tree_file += '%s & %s & %s &' % (num_ssms, num_cnvs, root.params)
  tree_file+='\\\\\n'

  for child in root.children():
    tree_file = print_index(child, tree_file)
  return tree_file

# writes the latex code
# tssb: tssb structure of the tree
# fout: output file for latex
def print_tree_latex(tssb,fout):
  global count
  #remove_empty_nodes(tssb.root, None)

  fout = open(fout,'w')
  count=-1
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
  tree_file=print_index(tssb.root['node'], tree_file)
  tree_file+='\\hline\n'
  tree_file+='\\end{tabular}\n'
  tree_file+='};\n'
  tree_file+='\\end{tikzpicture}\n'
  tree_file+='};\n'  
  tree_file+='\\node at (b.south) [anchor=north,yshift=-.5cm]{Posterior probability: 0};\n'
  tree_file+='\\end{tikzpicture}\n'
  tree_file+='\end{document}\n'
  fout.write(tree_file)
  fout.close()  
    
if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    description='Plot posterior trees resulting from PhyloWGS run',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('ssm_file',
    help='File listing SSMs (simple somatic mutations, i.e., single nucleotide variants. For proper format, see README.txt.')
  parser.add_argument('cnv_file',
    help='File listing CNVs (copy number variations). For proper format, see README.txt.')
  parser.add_argument('trees_dir',
    help='Directory where the MCMC trees/samples are saved')
  args = parser.parse_args()

  try:
    os.mkdir('latex')
  except OSError, e:
    if e.errno == 17: # Directory exists
      pass
    else:
      raise e

  compute_lineages(args.trees_dir, args.ssm_file, args.cnv_file)
