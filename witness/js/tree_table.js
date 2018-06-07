function TreeTable(){
}

TreeTable.prototype._predetermined_permutation = function(arr,n) {
  //This is used to choose n random trees from the unclustered trees. Because
  //I don't want the trees displayed to change every time we switch windows,
  //I will make use of the predetermined random numbers in 
  //Config.tree_viewer.tree_index_determinants.
  for(var idx = 0; idx < arr.length; idx++){
    var swpIdx = idx + Math.floor(Config.tree_viewer.tree_index_determinants[idx] * (arr.length - idx));
    // now swap elements at idx and swpIdx
    var tmp = arr[idx];
    arr[idx] = arr[swpIdx];
    arr[swpIdx] = tmp;
  }
  if(typeof n == "undefined" | n >= arr.length){
    return arr;
  }else{
    ans = [];
    for(var idx = 0; idx < n; idx++){
      ans.push(arr[idx]);
    }
    return ans;
  }
}

TreeTable.prototype._get_trees_from_clusters = function(trees, clusters){
  //Show a specified number of trees from each cluster. One of these should be the
  //representative tree for that cluster, and the others will be picked at random.
  var tree_indices = [];
  var rand_num_idx = 0;
  Object.keys(clusters).forEach(function(cidx){
    clust_members = clusters[cidx].members;
    if (clust_members.length<=Config.tree_viewer.num_trees_to_show){
      tree_indices = tree_indices.concat(clust_members);
    }else{
      //Representative tree should always be shown
      var rep_tree = clusters[cidx].representative_tree;
      tree_indices.push(rep_tree);
      //Get the indices of the other random trees. They should not repeat nor should they contain the rep_tree index
      var trees_added = 0;
      while(trees_added < Config.tree_viewer.num_trees_to_show-1){
        this_mem = Math.floor( Config.tree_viewer.tree_index_determinants[rand_num_idx] * clust_members.length);
        if(!(tree_indices.includes(clust_members[this_mem])) && !(clust_members[this_mem] == rep_tree)){
          tree_indices.push(clust_members[this_mem]);
          trees_added = trees_added + 1;
        }
        rand_num_idx = rand_num_idx + 1;
      }
    }
  })
  return tree_indices;
}
  
TreeTable.prototype._determine_clusters_to_report = function(trees, clusters, have_clust_data, report_small_clusters){
  if(!have_clust_data)
    return {to_report: clusters, rest: []};
  var separated_clusters = ClusterUtil.separate_clusters_by_size(clusters, Object.keys(trees).length*Config.small_cluster_tree_prop_cutoff);
  var reported_clusters = report_small_clusters ? clusters : separated_clusters.large;
  var unreported_clusters = report_small_clusters ? [] : separated_clusters.small;

  return {reported: reported_clusters, rest: unreported_clusters};
}
  
TreeTable.prototype._determine_trees_for_table = function(trees, reported_clusters, unreported_clusters, have_clust_data, show_all_trees){
  //If we don't have cluster data then it only makes sense to show all
  //trees and report that they don't belong to a cluster.
  if(!have_clust_data){
    var tree_indices = Util.sort_ints(Object.keys(trees).map(function (a) {return parseInt(a,10)}));
    return tree_indices;
  }
  // Get the indices of the trees and the clusters they belong to
  if(show_all_trees){
    var tree_indices = Util.sort_ints(Object.keys(trees).map(function (a) {return parseInt(a,10)}));
  }else{
    var clustered_trees = this._get_trees_from_clusters(trees, reported_clusters);
    var unclustered_trees = this._get_trees_from_clusters(trees, unreported_clusters);
    unclustered_trees = this._predetermined_permutation(unclustered_trees,Config.tree_viewer.num_trees_to_show);
    var tree_indices = (clustered_trees).concat(unclustered_trees);
  }
  return tree_indices;
}
  
TreeTable.prototype.create = function(summary, container) {
  //Based on whether or not we have cluster information, and on the settings in Config.tree_viewer,
  //create and populate the tree table in tree viewer with all required information.

  //Determine the factors that will determine what we show.
  var have_clust_data = Util.have_cluster_data(summary);
  var report_small_clusters = Config.report_small_clusters;
  var show_all_trees = Config.tree_viewer.show_all_trees;
  //Determine what clusters we are reporting, if any, and the trees that will populate the table.
  var clusters = this._determine_clusters_to_report(summary.trees, summary.clusters, have_clust_data, report_small_clusters);
  var tree_indices = this._determine_trees_for_table(summary.trees, clusters.reported, clusters.rest, have_clust_data, show_all_trees);
  //Roundabout way of calcing num_samples
  var first_tree_idx = tree_indices[0];
  var first_pop_idx = Object.keys(summary.trees[first_tree_idx].populations)[0];
  var num_samples = summary.trees[first_tree_idx].populations[first_pop_idx].cellular_prevalence.length;
  //Fill in the tree table.
  tree_indices.forEach(function(tidx) {
    //calc trees llh
    var total_ssms = 0;
    Object.keys(summary.trees[tidx].populations).forEach(function(pidx) {
      total_ssms += summary.trees[tidx].populations[pidx].num_ssms;
    })
    var normllh_nats = -summary.trees[tidx].llh / (total_ssms * num_samples);
    var normllh_bits = normllh_nats / Math.log(2);
    //Determine which cluster the tree belongs to, if any.
    var cluster = have_clust_data ? TreeUtil.find_cluster_from_treeidx(tidx, clusters.reported): -1;
    //Finally, fill in a row with the information we have about the tree.
    var row = '<td class="tree-index">' + tidx + '</td>'
      + '<td class="tree-llh">' + normllh_bits.toFixed(1) + '</td>'
      + '<td class="tree-nodes">' + Object.keys(summary.trees[tidx].populations).length + '</td>'
      + (cluster==-1 ? '<td class="cluster" data-sort-value="-1">-</td>' : '<td class="cluster">' + cluster + '</td>');
    ['linearity_index', 'branching_index', 'clustering_index'].forEach(function(idxname) {
      var val = summary.trees[tidx].hasOwnProperty(idxname) ? summary.trees[tidx][idxname].toFixed(2) : '&mdash;';
      row += '<td>' + val + '</td>';
    })
    $('<tr/>').html(row).appendTo(container);
  })
  return;
}