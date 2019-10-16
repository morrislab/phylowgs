import json
import pandas as pd
import os
from collections import Counter
import zipfile
import argparse

def return_cnv_table(cnv_data):
    full_df = cnv_data.drop(['a', 'd'], axis=1)
    one = cnv_data['physical_cnvs'].str.split(",",expand = True)
    one = one.replace({"[^0-9|.]": ''}, regex=True)
    full_df = full_df.join(one,how = "outer")
    full_df = full_df.rename(columns={ 0: "chr", 1: "start",2: "end", 3: "major_cn", 4 : "minor_cn", 5 : 'cell_prev','cnv':'id'})
    full_df = full_df.drop("physical_cnvs",axis = 1)
    return(full_df)

def return_ssm_table(ssm_data):
    
    full_df = ssm_data.drop(['mu_r', 'mu_v'], axis=1)
    pd.set_option('mode.chained_assignment', None)
    full_df['chr'], full_df['start'] = full_df['gene'].str.split('_', 1).str
    full_df = full_df.drop("gene", axis = 1)
   
    try:
        one = full_df['a'].str.split(",",expand = True)
        two = full_df['d'].str.split(",",expand = True)
        one.columns = [str(col) + '_ref' for col in one.columns]
        two.columns = [str(col) + '_total' for col in two.columns]
        full_df = full_df.join([one, two],how = "outer")    
    except:
        full_df = full_df.rename(columns={'a' : 'ref', 'd' : 'total'})
    #one.columns = [str(col) + '_ref' for col in one.columns]
    #two.columns = [str(col) + '_total' for col in two.columns]
    #full_df = full_df.join([one, two],how = "outer")
    #full_df = full_df.drop(['a', 'd'], axis = 1)
    return(full_df)

def return_top_k_trees(summary_file, k = 5):
    with open(summary_file, "r") as read_file:
        data_summary = json.load(read_file)
    trees = {}
    for tr in data_summary['trees']:
        trees[tr] = (data_summary['trees'][tr]['llh'])
    c = Counter(trees)
    
    return(c.most_common(k))

def tree_json(mutasgn_path, tree):
    mutzip = zipfile.ZipFile(mutasgn_path)

    tree_json = mutzip.open((tree + ".json"))
    contents = tree_json.read()
    tree_json.close()
    #parse the json so it's a sensible looking json
    my_json = contents.decode('utf8').replace("'", '"')
    data = json.loads(my_json)
    return(data)

def tree_to_df(tree_data, tree_number):
    #go through all the populations
    df = pd.DataFrame()
    for p in tree_data['mut_assignments']:
        df2 = pd.DataFrame(tree_data['mut_assignments'][p]['ssms'])
        df2['population'] = p
        df2['aberation'] = 'ssm'
        df3 = pd.DataFrame(tree_data['mut_assignments'][p]['cnvs'])
        df3['population'] = p
        df3['aberation'] = 'cnvs'
        df4 = pd.concat([df2,df3],sort=False)
        df = pd.concat([df,df4],sort=False)
    df['tree'] = tree_number

    df.columns = ['id','population','aberation','tree']

    return(df)
    
def read_best_trees_json(mutasgn_path, summary_file, k0 = 5):
    #get the tree with the max llh, as per Morris comment that tree with max llh is 'best'
    top_k_trees = return_top_k_trees(summary_file, k = k0)
    #use the zip file module to read in the zip file
    #get the json of the best tree
    best_trees = [x[0] for x in top_k_trees]
    all_trees = []
    for tree in best_trees:
        tree_data = tree_json(mutasgn_path, tree)
        tree_df = tree_to_df(tree_data, tree)
        all_trees.append(tree_df)
    return(all_trees)

def add_labels_best_trees(all_trees, all_data_simple):
    merge_trees = []
    for tree in all_trees:
        tree = pd.merge(tree, all_data_simple, on = 'id')
        merge_trees.append(tree)
    return(merge_trees)
    
def write_tree_likelihoods(output_folder, summary_file, k0 = 5):
    tree_likes = return_top_k_trees(summary_file, k = k0)
    with open((output_folder + "/" + 'tree_likelihoods.txt'), 'w') as f:
        f.write('\n'.join('%s %s' % x for x in tree_likes))

def write_tree_csvs(output_folder, all_trees):
    for tr in all_trees:
        tr.to_csv(output_folder + "/" + str((tr['tree'].iloc[0])) + ".csv")
    

def main():

    parser = argparse.ArgumentParser(
    description='Parse the output of PhyloWGS into useful tables for downstream analysis and exploration',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--cnv_input', dest='cnv_input', action='append', required=True,
    	help='The cnv.txt file which was originally fed into PhyloWGS')
    parser.add_argument('--ssm_input', dest='ssm_input', action='append', required=True,
    	help='The ssm.txt file which was originally fed into PhyloWGS')
    parser.add_argument('--summary_file', dest='summary_file', action='append', required=True,
    	help='File path to PhyloWGS output *.summ.json file, should be unziped with gunzip first')
    parser.add_argument('--mutasgn_path', dest='mutasgn_path', action='append', required=True,
    	help='File path to PhyloWGS output *.mutass.zip file')
    parser.add_argument('--output_folder', dest='output_folder', action='append', required=True,
    	help='Folder you would like the parsed output sent to')
    parser.add_argument('--k', dest='k', type=int, default=5,
    	help='K for how many of the top tress to output')

    args = parser.parse_args()
 
    cnv_simple = return_cnv_table(pd.read_csv(args.cnv_input[0], sep = "\t"))
    #print("read the cnv")
    #print(cnv_simple.head())
    
    #print(args.ssm_input[0])
    temp = args.ssm_input[0]
    ssm_simple = return_ssm_table(pd.read_csv(temp, sep = "\t"))
    all_data_simple = pd.concat([cnv_simple,ssm_simple],sort = False)

    all_trees = read_best_trees_json(args.mutasgn_path[0], args.summary_file[0], k0 = args.k)
    all_trees = add_labels_best_trees(all_trees, all_data_simple)

    write_tree_csvs(args.output_folder[0], all_trees)
    write_tree_likelihoods(args.output_folder[0], args.summary_file[0], k0 = args.k)

if __name__ == '__main__':
  main()
