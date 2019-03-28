import numpy as np
#import scipy as sp
from sklearn.metrics import jaccard_similarity_score
from IPython import embed
import time

class SSM_Analyser:
  """
    TODO:
    - Want to calculate three stats:
      1) DONE What is the cluster-concensus assignment and confidence for each SSM. Assignment
      to be determined by determining which cluster the SSM went into the most and the 
      confidence will be how many times it went into that cluster.
      2) What is the concensus phi for the cluster? Based on matching clusters across trees
      should get a distribution of phis. Can report the mean and condifence interval of 
      this distribution. Will have to weight each phi by the number of SSMs beloning to the
      cluster.
      3) What is the mean number of clusters from each tree that match to a given cluster in R?
      E.g., if this is around 2, it says that the reference tree merged clusters that are 
      usually split in other trees.
    - Write the results to a file, no need to integrate into witness right off the bat. Plus
    I probably don't have time anyways.
    - Match clusters between trees using the jaccard score.
    - Only compare trees with the representative tree for the cluster they are in.
    - probably rename this...
    NOTATION NOTE: I will call nodes the clusters of ssms within a tree, I will call a cluster a cluster of trees.
    """
  def analyse(self,clusters,summaries,mutlist,mutass):
    """
    :parameters: clusters, which is a dictionary containing information on clusters.
                 mutlist, which is a dictionary containing information on the ssms
                 mutass, which contains tells us which ssms belong to which clusters for each tree
    :results: for each cluster of tree, which ssms go to which clusters most often, what the concencus phi 
              if for each cluster, and the mean number of clusters from each tree that mathes with a given 
              cluster in R.
    """

    clusters = self._organize_trees_by_clusters(clusters, summaries, mutass, mutlist)

    ssm_stats = self._calc_ssm_stats(clusters, mutlist)

    #Calculate the concensus phi for each cluster in the rep tree using the phis from all other trees. 
    #Note: Ensure each phi is weighted by the number of ssms in the cluster
    consensus_phis = self._calc_consensus_phis(clusters, summaries)

    #Calc the mean number of clusters that maps to the reference cluster. Should tell us how often a cluster is split or merged with another.
    mean_num_clusts = self._calc_mean_num_node_mappings(clusters, summaries)

    self._write_results(ssm_stats, consensus_phis, mean_num_clusts)

  def _organize_trees_by_clusters(self, clusters, summaries, mutass, mutlist):
    org_clusters = {}
    for cidx in clusters:
      rep_tree_idx = clusters[cidx]['representative_tree']
      org_clusters[cidx] = {}
      org_clusters[cidx]['representative_tree'] = rep_tree_idx
      ref_tree_ssm_dist = mutass[rep_tree_idx]
      trees = {}
      for midx in clusters[cidx]['members']:
        trees[midx] = {
          "mut_assignments": mutass[midx],
          "num_nodes": len(summaries[midx]['populations']),
          "structure": summaries[midx]['structure'],
          "node_mapping_to_ref": self._match_clusters_to_rep_tree(ref_tree_ssm_dist, mutass[midx], mutlist['ssms'].keys())
        }
      org_clusters[cidx]['members'] = trees
    return org_clusters

  def _match_clusters_to_rep_tree(self, ref_tree_ssm_dist, comp_tree_ssm_dist, all_ssms):
    # first, for each node create a boolean array specifying whether or not an ssm is in a node
    #This step takes forever (9min for 2500 trees and 2000 ssms). See if can make it faster
    mapping = {};
    comp_nodes = comp_tree_ssm_dist.keys()
    ref_nodes  = ref_tree_ssm_dist.keys()
    for comp_node in comp_nodes:
      comp_ssms = comp_tree_ssm_dist[comp_node]['ssms']
      comp_bool = [(x in comp_ssms) for x in all_ssms]
      jacs = [];
      for ref_node in ref_nodes:
        ref_ssms = ref_tree_ssm_dist[ref_node]['ssms']
        ref_bool = [(x in ref_ssms) for x in all_ssms]
        jacs.append(jaccard_similarity_score(ref_bool, comp_bool))
      mapping[comp_node] = comp_nodes[np.argmax(jacs)]
    return mapping

  def _calc_ssm_stats(self, clusters, mutlist):
    #Create a dict of ssms, each saying which node they belong to in the ref tree, and what our confidence is on that.
    #SIGH, NEED TO DO THIS FOR EVERY CLUSTER INIVIDUALLY, OTHERWISE NODE NUMBER COULD REFER TO DIFFERENT NODES
    embed()
    
    ssm_info = {}
    ssm_info['cluster'] = {}
    for cidx in clusters.keys():
      cluster = clusters[cidx]
      nTrees = len(clusters[cidx]['members'])
      ssm_info['cluster'][cidx] = {}
      for ssm in mutlist['ssms'].keys():
        ssm_info['cluster'][cidx][ssm] = {}
        ssm_info['cluster'][cidx][ssm]['nodes'] = []
      for tree in cluster['members'].values():
        for node in tree['mut_assignments'].keys():
          node_ssm_assignments = tree['mut_assignments'][node]['ssms']
          for ssm in node_ssm_assignments:
            ssm_info['cluster'][cidx][ssm]['nodes'].append(tree['node_mapping_to_ref'][node])

      for ssm in mutlist['ssms'].keys():
        ssm_info['cluster'][cidx][ssm]['assignment'] = max(set(ssm_info['cluster'][cidx][ssm]['nodes']), key=ssm_info['cluster'][cidx][ssm]['nodes'].count)
        ssm_info['cluster'][cidx][ssm]['confidence'] = ssm_info['cluster'][cidx][ssm]['nodes'].count(ssm_info['cluster'][cidx][ssm]['assignment']) / float(nTrees)
    return ssm_info

  def _calc_consensus_phis(self, clusters, summaries):
    #this will return, for each reference tree in each cluster, the consensus phis for each node, calculated
    #by summing the phi for the nodes that map to it and dividing by the total number of trees in that cluster.
    embed()
    nSamples = len(summaries[0]['populations'][0]['cellular_prevalence'])
    consensus_phis = {}
    for cidx in clusters.keys():
      cluster = clusters[cidx]
      rep_tree = cluster['representative_tree']
      n_members = len(cluster['members'])
      consensus_phis[cidx] = {}
      for p in summaries[rep_tree]['populations'].keys():
        consensus_phis[cidx][p] = np.zeros(nSamples)
      for midx in cluster['members'].keys():
        tree_info = cluster['members'][midx]
        for this_node in tree_info['node_mapping_to_ref'].keys():
          ref_node = tree_info['node_mapping_to_ref'][this_node]
          consensus_phis[cidx][ref_node] += np.array(summaries[midx]['populations'][this_node]['cellular_prevalence'])
      for p in summaries[rep_tree]['populations'].keys():
        consensus_phis[cidx][p] = (consensus_phis[cidx][p]/float(n_members)).tolist()
    return consensus_phis

  def _calc_mean_num_node_mappings(self, clusters, summaries):
    #Go through the mappings and see how many nodes are, on average, mapped to a node in the reference tree
    embed()
    mean_num_clusts = [];
    for cluster in clusters.values():
      rep_tree = cluster["representative_tree"]
      n_rep_nodes = len(summaries[rep_tree]['populations'])
      num_mappings = np.zeros(n_rep_nodes)
      for tree in cluster['members'].values():
        for mapped_to in tree['node_mapping_to_ref'].values():
          num_mappings[mapped_to] += 1
      mean_num_clusts.append((num_mappings/len(cluster['members'])).tolist())
    return mean_num_clusts

  def _write_results(self, ssm_stats, consensus_phis, mean_num_clusts):
    
    
    return