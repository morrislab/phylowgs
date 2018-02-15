import numpy as np
from sklearn.mixture import GaussianMixture
from pwgsresults.index_calculator import IndexCalculator
from scipy import linalg

class GMMAnalyzer:
  def analyze(self,summaries):
    """
    Takes the tree summaries from ResultGenerator().generate, calculates the linear, branching and clustering indicies
    for them, and then clusters the trees based on their linearity and branching indicies.
    :param summaries: the tree summaries
    :return: formated dictionary containing all relevant information for the clusters, ready to be inserted to the .summ json
    """
    data = self._calc_gmm_input(summaries)
    gmm_run_info = self._run_gmm(data);
    clustInfo = self._generateOutput(gmm_run_info, summaries);

    return clustInfo
  
  
  def _calc_gmm_input(self,summs):
    """
    Makes use of IndexCalculator class to calculate the linear, branching and clustering indicies for use in gmm clustering.
    :param summs: the summaries of the trees as output by ResultGenerator().generate()
    :return: two sets of inputs for gmm_run, one where we are clustering using just LI, BI, and another using BI/(LI+BI), CI
    """
    out = {
      "LI_BI": [],  #[LI,BI]
      "CI_nBI": [] #[BI/(LI+BI),CI],
    }
    for summary in summs.values():
      calculator = IndexCalculator(summary)
      LI = summary['linearity_index'] = calculator.calc_linearity_index()
      BI = summary['branching_index'] = calculator.calc_branching_index()
      CI = summary['clustering_index'] = calculator.calc_clustering_index()
      out["LI_BI"].append([LI, BI])
      out["CI_nBI"].append([CI, BI/(LI+BI)])
    
    out["LI_BI"] = np.array(out["LI_BI"])
    out["CI_nBI"] = np.array(out["CI_nBI"])
    return out
  
  
  def _run_gmm(self,data):
    """
    Runs gmm using BIC to determine number of components, and return associated data
    :param data: LI/BI coordinates for all sampled trees
    :return: Weight, mean, covariance for each cluster, and assignments for each sampled tree
    """
    out = {};
    for key in data:
      thisdata = data[key]
      num_components = self._get_components_min_bic(thisdata)
      gmm = GaussianMixture(n_components=num_components, n_init=2, covariance_type="full").fit(thisdata)
      out[key] = {
        "weights": gmm.weights_,
        "means": gmm.means_,
        "covariances": gmm.covariances_,
        "members": gmm.predict(thisdata),
        "responsibilities": gmm.predict_proba(thisdata),
        "ellipses": self._generateEllipseInfo(gmm.means_,gmm.covariances_)
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

  def _generateEllipseInfo(self,means,covs):
    """
    Take the output from GMM analysis and generate a dictionary whose elements describe an ellipse 
    and can be used by witness to plot the clusters
    :params means, covs: The means and covariance matrix output from gmm analysis
    :return: dictionary containing ellipse mean, major_axis, minor_axis, and rotation
    """
    ellDict = [];
    for i in range(len(means)):
        # Note: v is the eigenvalues of the covs matrix. These two values represent the highest variance and the variance of the vector orthogonal to the vector with the highest variance (which is also the lowest variance)
        # Note: w is the eigenvectors of the covs matrix. W[0,:] is the vector representing the direction of highest variance, and W[1,:], the lowest variance.
        v, w = linalg.eigh(covs[i])
        # Take 2sqrt(2v) to convert v to represent 2 standard deviations from the mean
        v = np.multiply(2.,  np.sqrt(2.) * np.sqrt(v))
        angle = np.arctan(w[0,1] / w[0,0])
        ellDict.append({
            "mean": means[i].tolist(),
            "angle": angle,
            "major_axis":  v[0],
            "minor_axis": v[1]
        })
    return ellDict


  def _generateOutput(self, gmm_run_info, summaries):
    """
    Take the results of GMM analysis and create a dictionary that can be inserted into the .summ file
    :params weights, means, covs, assignments, responsibilities: results from gmm_run
    :return: organized dictionary ready for output to .summ json
    """
    out = {}
    for key in gmm_run_info:
      thisrun = gmm_run_info[key];
      out[key] = {};
      nclust = len(thisrun["weights"]);
      for i in range(nclust):
        out[key][str(i)] =  {
          "weight": thisrun["weights"][i],
          "members": [j for j,k in zip(summaries.keys(), thisrun["members"]) if k==i], 
          "responsibilities": thisrun["responsibilities"][:,i].tolist(),
          "mean": thisrun["means"][i].tolist(),
          "covariance": thisrun["covariances"][i].tolist(),
          "ellipse": thisrun["ellipses"][i]
          };
    return out
    
    
    
    
    
    
    
    
    
    
    

