import numpy as np
from sklearn.mixture import BayesianGaussianMixture, GaussianMixture
from pwgsresults.index_calculator import IndexCalculator
from pwgsresults.spectral_clustering import SpectralClustering
from scipy import linalg
import scipy as sp

class TreeClusterer:
  def find_clusters(self,summaries):
    """
    This clusters the trees based on their structure.
    First, this takes the tree summaries from ResultGenerator().generate() and calculates the linearity, branching and clustering indexes.
    Second, we cluster the linear trees (ie, branching index ==0) based on the number of nodes that they have.
    Third, we cluster the remaining trees based on the indexes calculated in the first step.
    Fourth, output the cluster information in a json-friendly format.
    Note that currently we clustering using two metrics: [LI, BI] and [CI, BI/(LI+BI)], so the final output dictionary has a different entry for each metric. 
    :param summaries: the tree summaries output from ResultGenerator().generate()
    :return: formated dictionary containing all relevant information for the clusters, ready to be inserted to the .summ json
    """
    self._calc_tree_structure_indexes(summaries)
    is_linear, clustering_data, tree_idx, n_nodes = self._get_tree_clustering_info(summaries)
    clustInfo = self._determine_clusters(tree_idx, clustering_data, n_nodes, is_linear)
    return clustInfo
  
  def _calc_tree_structure_indexes(self,summs):
    """
    Makes use of IndexCalculator class to calculate the linear, branching and clustering indexes and saves the indexes in the summary dictionary. Returns None
    :param summs: the summaries of the trees as output by ResultGenerator().generate()
    "return: none
    """
    for summary in summs.values():
      calculator = IndexCalculator(summary)
      summary['linearity_index'] = calculator.calc_linearity_index()
      summary['branching_index'] = calculator.calc_branching_index()
      summary['clustering_index'] = calculator.calc_clustering_index()
  
  def _get_tree_clustering_info(self,summs):
    """
    Takes the tree summaries and extracts the relevant tree information that will be used for clustering.
    :param summs:the summaries of the trees as output by ResultGenerator().generate()
    :return is_linear: list of booleans answering whether or not a tree at that index is linear or not
    :return clustering_data: list of data points to be used for clustering the trees. distance metric is [CI,BI/(LI+BI)]
    :return tree_idx: list of the indexes assigned to each tree.
    :return n_nodes: list representing the number of nodes found for each tree.
    """
    
    is_linear = [summ["branching_index"]==0 for summ in summs.values()]
    clustering_data = [[summ["clustering_index"], summ["branching_index"] / (summ["branching_index"] + summ["linearity_index"])] for summ in summs.values()]
    tree_idx = [tree_idx for tree_idx in summs.keys()]
    n_nodes = [len(summ["populations"]) for summ in summs.values()]
    
    return is_linear, clustering_data, tree_idx, n_nodes
  
  def _determine_clusters(self, tree_idxs, clustering_data, num_nodes, is_linear):
    """
    This will split the data into trees that are linear and not-linear and then run the proper analysis for each set.
    :param tree_idxs: a list of the indexes of the trees. 
    :param clustering_data: list of data points from each tree.
    :param num_nodes: list of he number of nodes associated with each tree
    :param is_linear: list of booleans answering whether or not a tree is linear.
    :return: dictionary which represents the clusters found in this anaysis. linear and non-linear trees are merged into the same dictionary.
    """
    
    #First, let's cluster the non-linear trees. 
    non_lin_data = [x for x,this_is_lin in zip(clustering_data,is_linear) if not this_is_lin]
    non_lin_tree_idxs = [x for x,this_is_lin in zip(tree_idxs,is_linear) if not this_is_lin]
    non_lin_clusters = self._run_spectral_clustering(non_lin_data, non_lin_tree_idxs)
    
    #Now cluster the linear trees based on the number of nodes they have
    lin_data = [x for x,this_is_lin in zip(clustering_data, is_linear) if this_is_lin]
    lin_tree_idxs = [x for x,this_is_lin in zip(tree_idxs, is_linear) if this_is_lin]
    lin_num_nodes = [x for x,this_is_lin in zip(num_nodes, is_linear) if this_is_lin]
    #For linear trees the minor axis for the ellipse is arbitrary as all y values are 0. So I choose a value here.
    ellipse_minor_axis = 1/10 * max([ datum[1] for datum in clustering_data ]) 
    if ellipse_minor_axis==0: ellipse_minor_axis = 0.01
    lin_clusters = self._calc_lin_tree_clusters(lin_data, lin_tree_idxs, lin_num_nodes, ellipse_minor_axis) 

    non_lin_clusters.update(lin_clusters)
    return non_lin_clusters
  
  def _run_spectral_clustering(self,data,tree_idxs):
    """
    Clusters the given tree index data using spectral clustering. See sklearn
    :param data: index coordinates for all sampled trees
    :param tree_idxs: the tree indexes corresponding to each tree input into data.
    :return: Weight, mean, covariance for each cluster, assignments for each sampled tree, and ellipse information that describes the tree and can be used for plotting.
    """
    #If all trees are linear then the inputs will be empty. Return an empty dictionary
    if not data:
      return {}

    data = np.array(data)
    max_num_clusters = 6 
    spectral = SpectralClustering(n_clusters = range(2,max_num_clusters+1)).fit(data)
    #Create an output dictionary that will contain all of the relevant information for each cluster.
    out = {}
    cluster_assignments = spectral.labels_
    num_clusters = spectral.n_clusters
    for this_clust_idx in range(0,num_clusters):
      this_clust_member_idxs = [idx for idx, cluster_assignment in zip(range(len(tree_idxs)), cluster_assignments) if cluster_assignment==this_clust_idx]
      this_clust_members_tree_idxs = [tree_idx for tree_idx, cluster_assignment in zip(tree_idxs, cluster_assignments) if cluster_assignment==this_clust_idx]
      this_clust_rep_idx = self._determine_representative_tree(data[this_clust_member_idxs,0],data[this_clust_member_idxs,1])
      rep_tree_idx = this_clust_members_tree_idxs[this_clust_rep_idx]
      out[str(this_clust_idx)] =  {
        "is_linear": False,
        "members": this_clust_members_tree_idxs,
        "representative_tree": rep_tree_idx
        }
    return out
  
  def _calc_lin_tree_clusters(self, data, tree_idxs, num_nodes, ellipse_minor_axis):
    """
    Clusters the linear trees based on the number of nodes that they have and returns a dictionary in the same format as run_gmm output
    :param data: list of data points from each tree.
    :param tree_idxs: a list of the indexes of the trees. 
    :param num_nodes: list of the number of nodes associated with each tree
    :param ellipse_minor_axis: minor axis of the ellipse that will describe the tree. Can't find this from just the linear trees as they all lie on the line y=0. Instead, calc it and input it here.
    :return: dictionary in the same format as run_gmm output
    """
    #If there are no linear trees, return an empty dictionary
    if not data:
      return {}
    data = np.array(data)
    out = {}
    
    unique_num_nodes = list(set(num_nodes))
    for this_num_nodes in unique_num_nodes:
      #The indexes wrt all of the sampled trees, both linear and non-linear
      cluster_members = [tree_idx for tree_idx, this_tree_num_nodes in zip(tree_idxs, num_nodes) if this_tree_num_nodes==this_num_nodes]
      #The indexes wrt to the linear trees, inserted into this function
      this_data_member_idxs = [tree_idx for tree_idx, this_tree_num_nodes in zip(range(len(tree_idxs)), num_nodes) if this_tree_num_nodes==this_num_nodes]
      clust_rep_tree_idx = self._determine_representative_tree(data[this_data_member_idxs,0],data[this_data_member_idxs,1])
      rep_tree_idx = cluster_members[clust_rep_tree_idx]
      
      #Format the results to match the gmm results. For any fields that aren't necessary, set to None.
      out["linear_" + str(this_num_nodes)] =  {
        "is_linear": True,
        "members": cluster_members,
        "representative_tree": rep_tree_idx
        }
    return out

  def _determine_representative_tree(self,x,y):
    """
    Determine the best tree to represent the given members input (typically all tree members that belong to a cluster)
    """
    
    #Testing something. Because of the post-processing of the trees (deleting subclones) sometimes all trees in a cluster
    #will have the exact same branching, linearity and coclustering indicies. I will check for that and just say that, 
    #arbitrarily, the rep tree is the first in the cluster. They are all the same anyways so should all be "representative".
    #Still, should bring this up to Jeff.
    if all(i == x[0] for i in x) and all(i == y[0] for i in y):
      return 1

    data = np.vstack((x,y))
    try:
      density = list(sp.stats.gaussian_kde(data)(data))
    except (np.linalg.linalg.LinAlgError, FloatingPointError):
      density = list(sp.stats.gaussian_kde(x)(y))
    best_tree_idx = np.argmax(density)
    return best_tree_idx