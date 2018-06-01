function TreeViewer() {
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

TreeViewer.prototype._determine_table_trees = function(summary){
  var trees = summary.trees;
  var clusters = summary.clusters;
  if (Config.show_all_trees | !Util.have_cluster_data(summary)){
    return Util.sort_ints(Object.keys(trees).map(function (a) {return parseInt(a,10)}));
  }

  //Show a specified number of trees from each cluster. One of these should be the
  //representative tree for that cluster, and the others will be picked at random.
  var tree_indices = [];
  var rand_num_idx = 0;
  Object.keys(clusters).forEach(function(cidx){
    clust_members = clusters[cidx].members;
    if (clust_members.length<=Config.num_trees_to_show){
      tree_indices = tree_indices.concat(clust_members)
    }
    else{
      //Representative tree should always be shown
      rep_tree = clusters[cidx].representative_tree;
      tree_indices.push(rep_tree);
      //Get the indicies of 4 other random trees. They should not repeat nor should they contain the rep_tree index
      trees_added = 0;
      while(trees_added < Config.num_trees_to_show-1){
        this_mem = Math.floor( Config.tree_index_determinants[rand_num_idx] * clust_members.length );
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
  var tplotter = this;
  d3.json(dataset.summary_path, function(summary) {
    var have_clust_info = Util.have_cluster_data(summary);
    var separated_clusters = have_clust_info ? ClusterUtil.separate_clusters_by_size(summary.clusters, Object.keys(summary.trees).length*Config.tiny_cluster_criteria): -1;
    var tree_indices = tplotter._determine_table_trees(summary);
    tree_container.empty();
    var first_tree_idx = tree_indices[0];
    var first_pop_idx = Object.keys(summary.trees[first_tree_idx].populations)[0];
    var num_samples = summary.trees[first_tree_idx].populations[first_pop_idx].cellular_prevalence.length;

    tree_indices.forEach(function(tidx) {
      var total_ssms = 0;
      Object.keys(summary.trees[tidx].populations).forEach(function(pidx) {
        total_ssms += summary.trees[tidx].populations[pidx].num_ssms;
      });

      var normllh_nats = -summary.trees[tidx].llh / total_ssms;
      normllh_nats /= num_samples;
      var normllh_bits = normllh_nats / Math.log(2);
      var cluster = have_clust_info ? TreeUtil.find_cluster_from_treeidx(tidx, separated_clusters.large): "-";
      var row = '<td class="tree-index">' + tidx + '</td>'
        + '<td class="tree-llh">' + normllh_bits.toFixed(1) + '</td>'
        + '<td class="tree-nodes">' + Object.keys(summary.trees[tidx].populations).length + '</td>'
        + (cluster==-1 ? '<td class="cluster" data-sort-value="-1">None</td>' : '<td class="cluster">' + cluster + '</td>');
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
