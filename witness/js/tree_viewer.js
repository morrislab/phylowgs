function TreeViewer() {
}

TreeViewer.prototype._predetermined_permutation = function(arr,n) {
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

TreeViewer.prototype._plot_pop_vafs = function(dataset, tidx) {
  if(!(dataset.hasOwnProperty('muts_path') && dataset.hasOwnProperty('mutass_path')))
    return;

  var pop_vaf_plotter = new PopVafPlotter();
  d3.json(dataset.muts_path, function(muts) {
    d3.json(dataset.mutass_path + '/' + tidx + '.json', function(mutass) {
      pop_vaf_plotter.plot(muts, mutass);
    });
  });
}

TreeViewer.prototype._determine_table_trees = function(trees, clusters){
  //Show a specified number of trees from each cluster. One of these should be the
  //representative tree for that cluster, and the others will be picked at random.
  var tree_indices = [];
  var rand_num_idx = 0;
  Object.keys(clusters).forEach(function(cidx){
    clust_members = clusters[cidx].members;
    if (clust_members.length<=Config.tree_viewer.num_trees_to_show){
      tree_indices = tree_indices.concat(clust_members)
    }
    else{
      //Representative tree should always be shown
      var rep_tree = clusters[cidx].representative_tree;
      tree_indices.push(rep_tree);
      //Get the indicies of 4 other random trees. They should not repeat nor should they contain the rep_tree index
      var trees_added = 0;
      while(trees_added < Config.tree_viewer.num_trees_to_show-1){
        this_mem = Math.floor( Config.tree_viewer.tree_index_determinants[rand_num_idx] * clust_members.length );
        if(!(tree_indices.includes(clust_members[this_mem])) && !(clust_members[this_mem] == rep_tree)){
          tree_indices.push(clust_members[this_mem]);
          trees_added = trees_added + 1;
        }
        rand_num_idx = rand_num_idx + 1;
      }
    }
  })
  return tree_indices
}

TreeViewer.prototype.render = function(dataset) {
  $('#tree-list').show();
  var tree_container = $('#trees tbody');
  tree_container.empty();
  var tplotter = this;
  d3.json(dataset.summary_path, function(summary) {
    var have_clust_data = Util.have_cluster_data(summary);
    var show_all_trees = Config.tree_viewer.show_all_trees;
    var report_small_clusters = Config.report_small_clusters;
    if(have_clust_data){
        if(report_small_clusters){
            var clusters = summary.clusters;
            var tree_indices = show_all_trees | !have_clust_data
                ? Util.sort_ints(Object.keys(summary.trees).map(function (a) {return parseInt(a,10)}))
                : tplotter._determine_table_trees(summary.trees, clusters);
        }else{
            var separated_clusters = ClusterUtil.separate_clusters_by_size(summary.clusters, Object.keys(summary.trees).length*Config.small_cluster_tree_prop_cutoff);
            var non_clusters = separated_clusters.small;
            var clusters = separated_clusters.large;
            if(show_all_trees){
                var tree_indices = Util.sort_ints(Object.keys(summary.trees).map(function (a) {return parseInt(a,10)}));
            }else{
                var clustered_tree_indicies = tplotter._determine_table_trees(summary.trees, clusters);
                var unclustered_tree_indicies = tplotter._determine_table_trees(summary.trees, non_clusters);
                unclustered_tree_indicies = tplotter._predetermined_permutation(unclustered_tree_indicies,Config.tree_viewer.num_trees_to_show);
                var tree_indices = clustered_tree_indicies.concat(unclustered_tree_indicies);
            }
        }
    }else{
        var tree_indices = Util.sort_ints(Object.keys(summary.trees).map(function (a) {return parseInt(a,10)}));
    };
    var first_tree_idx = tree_indices[0];
    var first_pop_idx = Object.keys(summary.trees[first_tree_idx].populations)[0];
    var num_samples = summary.trees[first_tree_idx].populations[first_pop_idx].cellular_prevalence.length;
    // Fill in the tree table.
    tree_indices.forEach(function(tidx) {
      var total_ssms = 0;
      Object.keys(summary.trees[tidx].populations).forEach(function(pidx) {
        total_ssms += summary.trees[tidx].populations[pidx].num_ssms;
      });

      var normllh_nats = -summary.trees[tidx].llh / total_ssms;
      normllh_nats /= num_samples;
      var normllh_bits = normllh_nats / Math.log(2);
      var cluster = have_clust_data ? TreeUtil.find_cluster_from_treeidx(tidx, clusters): -1;
      var row = '<td class="tree-index">' + tidx + '</td>'
        + '<td class="tree-llh">' + normllh_bits.toFixed(1) + '</td>'
        + '<td class="tree-nodes">' + Object.keys(summary.trees[tidx].populations).length + '</td>'
        + (cluster==-1 ? '<td class="cluster" data-sort-value="-1">-</td>' : '<td class="cluster">' + cluster + '</td>');
      ['linearity_index', 'branching_index', 'clustering_index'].forEach(function(idxname) {
        var val = summary.trees[tidx].hasOwnProperty(idxname) ? summary.trees[tidx][idxname].toFixed(2) : '&mdash;';
        row += '<td>' + val + '</td>';
      });
      $('<tr/>').html(row).appendTo(tree_container);
    });

    $('#trees').stupidtable();

    var already_autosorted = false;
    $('#trees').bind('aftertablesort', function() {
      if(already_autosorted)
        return;
      // If any restore events are in progress (e.g., waiting on a timer),
      // don't automatically click on the first element.
      if(StateManager.restoring > 0)
        return;
      tree_container.find('tr:first').click();
      already_autosorted = true;
    });

    // If direction not specified, this can end up being ascending or
    // descending sort, depending on prior sort state of table.
    $('#tree-llh').stupidsort('asc');

    tree_container.find('tr').click(function(evt) {
      evt.preventDefault();
      var self = $(this);
      self.siblings().removeClass('active');
      self.addClass('active');

      var tidx = self.find('.tree-index').text();
      StateManager.update('tidx', tidx);
      var tree_plotter = new TreePlotter();

      tree_plotter.draw(summary.trees[tidx].populations, summary.trees[tidx].structure, summary.trees[tidx].root, summary.params);
      tplotter._plot_pop_vafs(dataset, tidx);
    });
    $('#tree-list').scrollTop(0);
  });
}
