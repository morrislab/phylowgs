import numpy as np
from sklearn.mixture import GaussianMixture
from pwgsresults.index_calculator import IndexCalculator
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
    self._calc_tree_structure_indexes(summaries);
    is_linear, LI_BI, CI_nBI, tree_idx, n_nodes = self._get_tree_clustering_info(summaries);
    clustInfo = {
        "LI_BI": self._determine_clusters(tree_idx, LI_BI, n_nodes, is_linear),
        "CI_nBI": self._determine_clusters(tree_idx, CI_nBI, n_nodes, is_linear)
    };
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
    :return LI_BI: list of data points to be used for the [LI,BI] distance metric.
    :return CI_nBI: list of data points to be used for the [CI,(BI/(LI+BI)] distance metric.
    :return tree_idx: list of indexes assigned to each tree
    :return n_nodes: list representing the number of nodes found for each tree.
    """
    
    is_linear = [summ["branching_index"]==0 for summ in summs.values()];
    LI_BI = [[summ["linearity_index"], summ["branching_index"]] for summ in summs.values()];
    CI_nBI = [[summ["clustering_index"], summ["branching_index"] / (summ["branching_index"] + summ["linearity_index"])] for summ in summs.values()];
    tree_idx = [tree_idx for tree_idx in summs.keys()];
    n_nodes = [len(summ["populations"]) for summ in summs.values()];
    
    return is_linear, LI_BI, CI_nBI, tree_idx, n_nodes
  
  def _determine_clusters(self, tree_idxs, data_to_cluster, num_nodes, is_linear):
    """
    This will split the data into trees that are linear and not-linear and then run the proper analysis for each set.
    :param tree_idxs: a list of the indexes of the trees. 
    :param data_to_cluster: list of data points from each tree.
    :param num_nodes: list of he number of nodes associated with each tree
    :param is_linear: list of booleans answering whether or not a tree is linear.
    :return: dictionary which represents the clusters found in this anaysis. linear and non-linear trees are merged into the same dictionary.
    """
    
    #First, let's cluster the non-linear trees. 
    non_lin_data = [x for x,this_is_lin in zip(data_to_cluster,is_linear) if not this_is_lin];
    non_lin_tree_idxs = [x for x,this_is_lin in zip(tree_idxs,is_linear) if not this_is_lin];
    non_lin_clusters = self._run_gmm(non_lin_data, non_lin_tree_idxs);
    
    #Now cluster the linear trees based on the number of nodes they have
    lin_data = [x for x,this_is_lin in zip(data_to_cluster, is_linear) if this_is_lin];
    lin_tree_idxs = [x for x,this_is_lin in zip(tree_idxs, is_linear) if this_is_lin];
    lin_num_nodes = [x for x,this_is_lin in zip(num_nodes, is_linear) if this_is_lin];
    #For linear trees the minor axis for the ellipse is arbitrary as all y values are 0. So I choose a value here.
    ellipse_minor_axis = 1/10 * max([ datum[1] for datum in data_to_cluster ]) 
    if ellipse_minor_axis==0: ellipse_minor_axis = 0.01;
    lin_clusters = self._calc_lin_tree_clusters(lin_data, lin_tree_idxs, lin_num_nodes, ellipse_minor_axis); 

    non_lin_clusters.update(lin_clusters)
    return non_lin_clusters
    
  def _run_gmm(self,data,tree_idxs):
    """
    Runs gmm using BIC to determine number of components, and return associated data
    :param data: index coordinates for all sampled trees
    :param tree_idxs: the tree indexes corresponding to each tree input into data.
    :return: Weight, mean, covariance for each cluster, assignments for each sampled tree, and ellipse information that describes the tree and can be used for plotting.
    """
    #If all trees are linear then the inputs will be empty. Return an empty dictionary
    if not data:
      return {};
    
    out = {};
    data = np.array(data);
    num_components = self._get_components_min_bic(data)
    gmm = GaussianMixture(n_components=num_components, n_init=2, covariance_type="full").fit(data)
    
    #Create an output dictionary that will contain all of the relevant information for each cluster
    cluster_assignments = gmm.predict(data)
    cluster_responsibilities = gmm.predict_proba(data)
    num_clusters = len(gmm.weights_);
    for clust_idx in range(num_clusters):
      clust_members = [tree_idx for tree_idx, cluster_assignment in zip(tree_idxs, cluster_assignments) if cluster_assignment==clust_idx];
      clust_rep_tree_idx = self._determine_representative_tree(data[clust_members,0],data[clust_members,1]);
      rep_tree_idx = clust_members[clust_rep_tree_idx];
      out[str(clust_idx)] =  {
        "num_nodes": None,
        "weight": gmm.weights_[clust_idx],
        "members": clust_members, 
        "responsibilities": cluster_responsibilities[:,clust_idx].tolist(),
        "mean": gmm.means_[clust_idx].tolist(),
        "covariance": gmm.covariances_[clust_idx].tolist(),
        "ellipse": self._generate_EMM_ellipse_info(gmm.means_[clust_idx], gmm.covariances_[clust_idx]),
        "representative_tree": rep_tree_idx
        };
    return out
  
  def _get_components_min_bic(self,data, end_early=False, delta=2):
    """
    Get the number of number of gmm components which result in lowest bic score
    :param data: LI/BI coordinates of all sampled trees
    :param end_early: End on asymptotic convergence (default False)
    :param delta: Percent change threshold for convergence
    :return: Number of components
    """
    min_bic = 0
    prev_bic = 10000
    bic_clusters = 1
    size = 50 if data.size/2 > 50 else int(data.size/2 - 1)

    for i in range(size):
      gmm = GaussianMixture(n_components=i+1, n_init=2, covariance_type='full').fit(data)
      bic = gmm.bic(data)
      # Check for convergence
      if end_early and (prev_bic/bic) - 1 > - delta:
        return i + 1
      elif bic < min_bic:
        bic_clusters = i+1
        min_bic = bic
      prev_bic = bic
    return bic_clusters
  
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
      return {};
    data = np.array(data);
    out = {};
    
    unique_num_nodes = list(set(num_nodes));
    for this_num_nodes in unique_num_nodes:
      #The indexes wrt all of the sampled trees, both linear and non-linear
      cluster_members = [tree_idx for tree_idx, this_tree_num_nodes in zip(tree_idxs, num_nodes) if this_tree_num_nodes==this_num_nodes];
      #The indexes wrt to the linear trees, inserted into this function
      this_data_member_idxs = [tree_idx for tree_idx, this_tree_num_nodes in zip(range(len(tree_idxs)), num_nodes) if this_tree_num_nodes==this_num_nodes];
      clust_rep_tree_idx = self._determine_representative_tree(data[this_data_member_idxs,0],data[this_data_member_idxs,1]);
      rep_tree_idx = cluster_members[clust_rep_tree_idx];
      cluster_responsibilities = [1*(this_tree_num_nodes == this_num_nodes) for this_tree_num_nodes in num_nodes]; #0 if this tree doesn't have the number of nodes, and 1s if it does
      cluster_mean = [np.mean(np.array([x[0] for x,nNodes in zip(data, num_nodes) if nNodes==this_num_nodes])), 0];
      xVals = [x[0] for x,n_nodes in zip(data, num_nodes) if n_nodes == this_num_nodes]
      yVals = [x[1] for x,n_nodes in zip(data, num_nodes) if n_nodes == this_num_nodes]
      
      #Format the results to match the gmm results. For any fields that aren't necessary, set to None.
      out["linear_" + str(this_num_nodes)] =  {
        "num_nodes": this_num_nodes,
        "weight": None,
        "members": cluster_members, 
        "responsibilities": cluster_responsibilities,
        "mean": cluster_mean,
        "covariance": None,
        "ellipse": self._generate_linear_ellipse_info(xVals,yVals,ellipse_minor_axis),
        "representative_tree": rep_tree_idx
        };
    return out
  
  def _generate_linear_ellipse_info(self, xVals, yVals, ellipse_minor_axis):
    
    ellipse_mean = [np.mean(xVals),0]
    ellipse_major_axis = max( [abs(max(xVals)-ellipse_mean[0]), abs(min(xVals)-ellipse_mean[0])] );
    
    ellDict = {
        "mean": ellipse_mean,
        "angle": 0,
        "major_axis": ellipse_major_axis,
        "minor_axis": ellipse_minor_axis
    }
    
    return ellDict
  
  def _generate_EMM_ellipse_info(self,mean,cov):
    """
    Take the output from GMM analysis and generate a dictionary whose elements describe an ellipse 
    and can be used by witness to plot the clusters
    :params means, covs: The means and covariance matrix output from gmm analysis
    :return: dictionary containing ellipse mean, major_axis, minor_axis, and rotation
    """
    # Note: v is the eigenvalues of the covs matrix. These two values represent the highest variance and the variance of the vector orthogonal to the vector with the highest variance (which is also the lowest variance)
    # Note: w is the eigenvectors of the covs matrix. W[0,:] is the vector representing the direction of highest variance, and W[1,:], the lowest variance.
    v, w = linalg.eigh(cov)
    # Take 2sqrt(2v) to convert v to represent 2 standard deviations from the mean
    v = np.multiply(2.,  np.sqrt(2.) * np.sqrt(v))
    angle = np.arctan(w[0,1] / w[0,0])
    ellDict = {
        "mean": mean.tolist(),
        "angle": angle,
        "major_axis":  v[0],
        "minor_axis": v[1]
    }
    return ellDict

  def _determine_representative_tree(self,x,y):
    """
    Determine the best tree to represent the given members input (typically all tree members that belong to a cluster)
    """
    data = np.vstack((x,y))
    try:
      density = list(sp.stats.gaussian_kde(data)(data));
    except (np.linalg.linalg.LinAlgError, FloatingPointError):
      density = list(sp.stats.gaussian_kde(x)(y));
    best_tree_idx = np.argmax(density);
    return best_tree_idx;