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
    indicies = self._calc_indicies(summaries)
    weights, means, covs, assignments, responsibilities = self._run_gmm(indicies);
    ellipseInfo = self._generateEllipseInfo(means,covs)
    clustInfo = self._generateOutput(weights, means, covs, assignments, responsibilities, summaries, ellipseInfo);
    return clustInfo
  
  
  def _calc_indicies(self,summs):
    """
    Makes use of IndexCalculator class to calculate the linear, branching and clustering indicies for use in gmm clustering.
    :param summs: the summaries of the trees as output by ResultGenerator().generate()
    """
    out = [];
    for summary in summs.values():
      calculator = IndexCalculator(summary)
      summary['linearity_index'] = calculator.calc_linearity_index()
      summary['branching_index'] = calculator.calc_branching_index()
      summary['clustering_index'] = calculator.calc_clustering_index() 
      out.append([summary[x] for x in ['linearity_index','branching_index']])
    return np.array(out)
  
  
  def _run_gmm(self,data):
    """
    Runs gmm using BIC to determine number of components, and return associated data
    :param data: LI/BI coordinates for all sampled trees
    :return: Weight, mean, covariance for each cluster, and assignments for each sampled tree
    """
    num_components = self._get_components_min_bic(data)
    gmm = GaussianMixture(n_components=num_components, n_init=2, covariance_type="full").fit(data)
    return gmm.weights_, gmm.means_, gmm.covariances_, gmm.predict(data), gmm.predict_proba(data)
  
  
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


  def _generateOutput(self, weights, means, covs, assignments, responsibilities, summaries, ellipseInfo):
    """
    Take the results of GMM analysis and create a dictionary that can be inserted into the .summ file
    :params weights, means, covs, assignments, responsibilities: results from gmm_run
    :return: organized dictionary ready for output to .summ json
    """
    nclust = len(weights);
    clusters = {}
    for i in range(nclust):
      clusters[str(i)] =  {
        "weight": weights[i],
        "members": [j for j,k in zip(summaries.keys(), assignments) if k==i], 
        "responsibilities": responsibilities[:,i].tolist(),
        "mean": means[i].tolist(),
        "covariance": covs[i].tolist(),
        "ellipse": ellipseInfo[i]
        };
    return clusters
    
    
    
    
    
    
    
    
    
    
    

