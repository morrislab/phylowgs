function ClusterViewer(){
}

ClusterViewer.prototype.render = function(dataset){
  self = this;
  d3.json(dataset.summary_path, function(summary) {
    if(Util.have_cluster_data(summary)){
      // These functions require that the summ file be created with the latest version of write_results.
      // If not, report in container that there is no cluster data and that write_results should be run again
      self._render_lin_idx_vs_branch_idx(summary);
      self._render_cluster_table(summary);
    }else{
      $('#container').text('No cluster data available. Rerun write_results.py to view tree clustering.');
    }
  })
}

ClusterViewer.prototype._render_lin_idx_vs_branch_idx = function(tree_summary) {
  if(!tree_summary.hasOwnProperty('tree_densities') | !tree_summary.hasOwnProperty('clusters')) {
    return;
  }

  // Calculcate the xData and yData from the indicies
  var clust_data = TreeUtil.calc_clustering_data(tree_summary.trees);
  var best_tree_idx = TreeUtil.find_best_tree(tree_summary.tree_densities);
  // Create the traces for the data points for clustered trees
  if(Config.report_small_clusters){
    var cluster_colours = this._determine_cluster_colour_ordering(tree_summary.clusters, Config.tree_summ.cluster_colours);
    var traces = this._create_cluster_scatter_traces(tree_summary.clusters, cluster_colours, clust_data.CI, clust_data.nBI);
  }else{
    var separated_clusters = ClusterUtil.separate_clusters_by_size(tree_summary.clusters, Object.keys(tree_summary.trees).length*Config.small_cluster_tree_prop_cutoff);
    var cluster_colours = this._determine_cluster_colour_ordering(separated_clusters.large, Config.tree_summ.cluster_colours);
    // Create the trace for the data points
    var traces = this._create_cluster_scatter_traces(separated_clusters.large, cluster_colours, clust_data.CI, clust_data.nBI);
    traces.push(this._create_tiny_cluster_scatter_trace(separated_clusters.small, clust_data.CI, clust_data.nBI));
  }
  //Use the traces and plotting options to actually plot our data.
  var layout = {
    title: "Clustering Degree vs. Branching Degree (best tree: " + best_tree_idx + ")",
    height: 1000,
    xaxis: { title: 'CI'},
    yaxis: { title: 'BI/(LI+BI)'},
    hovermode: 'closest',
    plot_bgcolor: '#440154',
    showlegend: true,
    legend: {x: 1.01,
             bgcolor: '440154',
             font: {size: 13, 
                    color: 'white'}}
  };
  var container = document.querySelector('#container');
  var plot_container = document.createElement('div');
  container.appendChild(plot_container);
  Plotly.newPlot(plot_container, traces, layout);
}
  
ClusterViewer.prototype._determine_cluster_colour_ordering = function(clusters,colors){
  // Will simply assign colours to clusters ordered from the top to the bottom.
  var self = this;
  var nCol = colors.length;
  var nClust = Object.keys(clusters).length;
  if (nCol>nClust){ // Since there are more colours than clusters, there is no need to order the colouring
    return colors.slice(0,nClust);
  }
  
  var clust_ypos = [];
  var clust_nums = [];
  var count = 0;
  Object.keys(clusters).forEach(function(cidx){
    clust_ypos.push(clusters[cidx].mean[1]);
    clust_nums.push(count);
    count = count + 1;
  })
  var ordered_colours = Array(nClust);
  for(count=0; count<nClust; count++){
    var cidx = clust_ypos.indexOf(Math.max(...clust_ypos));
    ordered_colours[clust_nums[cidx]] = colors[count % nCol];
    clust_ypos.splice(cidx,1)
    clust_nums.splice(cidx,1)
  }
  return ordered_colours;
}

ClusterViewer.prototype._create_tiny_cluster_scatter_trace = function(clusters, xData, yData) {
  var this_xData = [];
  var this_yData = [];
  var labels = [];
  var marker_symbols = [];
  var marker_sizes = [];
  var marker_colours = [];
  Object.keys(clusters).forEach(function(cidx){
    var members = clusters[cidx].members;
    Object.keys(members).forEach(function(midx){
      var tidx = members[midx];
      labels.push('Tree ' + tidx + ' - Unclustered');
      marker_symbols.push('dot');
      marker_sizes.push(6);
      marker_colours.push('white')
      this_xData.push(xData[tidx])
      this_yData.push(yData[tidx])
    });
  })
  var trace = {
    x: this_xData,
    y: this_yData,
    name: 'Unclustered',
    mode: 'markers',
    type: 'scatter',
    text: labels,
    showlegend: true,
    marker: {
      symbol: marker_symbols,
      size: marker_sizes,
      line: { width: 0 },
      colorscale: 'Viridis',
      color: marker_colours,
    }
  };
  return trace
}

ClusterViewer.prototype._create_cluster_scatter_traces = function(clusters, cluster_colours, xData, yData) {
  var traces = [];
  var colour_idx = 0;
  var label_clust_index = 1;
  Object.keys(clusters).forEach(function(cidx){
    var members = clusters[cidx].members;
    var this_xData = [];
    var this_yData = [];
    var labels = [];
    var marker_symbols = [];
    var marker_sizes = [];
    var marker_colours = [];
    var rep_tree_index = clusters[cidx].representative_tree;
    Object.keys(members).forEach(function(midx){
      var tidx = members[midx];
      labels.push('Tree ' + tidx + ' - Cluster ' + (label_clust_index));
      marker_symbols.push('dot');
      marker_sizes.push(6);
      if(rep_tree_index === tidx){
        marker_symbols.push('cross');
        marker_sizes.push(30);
      }
      marker_colours.push(cluster_colours[colour_idx])
      this_xData.push(xData[tidx])
      this_yData.push(yData[tidx])
    });
    colour_idx = colour_idx+1;
    traces.push({
      x: this_xData,
      y: this_yData,
      name: 'Cluster ' + label_clust_index,
      mode: 'markers',
      type: 'scatter',
      text: labels,
      showlegend: true,
      marker: {
        symbol: marker_symbols,
        size: marker_sizes,
        line: { width: 0 },
        colorscale: 'Viridis',
        color: marker_colours,
      }
    })
    label_clust_index = label_clust_index + 1;
  })
  return traces
}

ClusterViewer.prototype._render_cluster_table = function(summary) {
  var cluster_table = $('#cluster-table .tree-summary').clone().appendTo('#container');
  cluster_table = cluster_table.find('tbody');
  var separated_clusters = ClusterUtil.separate_clusters_by_size(summary.clusters, Object.keys(summary.trees).length * Config.small_cluster_tree_prop_cutoff)
  var clusters = Config.report_small_clusters ? summary.clusters : separated_clusters.large;
  var self = this;
  var clust_num = 0;
  Object.keys(clusters).forEach(function(cluster_idx) {
    clust_num = clust_num+1;
    var C = summary.clusters[cluster_idx];
    var propTrees = (C.members.length/Object.keys(summary.trees).length).toFixed(2);
    // Create a new container to hold the representative tree and then make the tree object. After the row is made, draw the tree in the container.
    var rep_tree_container = document.createElement('div');
    rep_tree_container.id = cluster_idx;
    var rep_tree_container_id = 'rep_tree_container-' + cluster_idx;
    var representative_tree_idx = C.representative_tree;

    var entries = [clust_num].concat(propTrees).concat(representative_tree_idx).concat('<div id='+rep_tree_container_id+'></div>').map(function(entry) {
      return '<td>' + entry + '</td>';
    }); 
    $('<tr/>').html(entries.join('')).appendTo(cluster_table);
    self._draw_rep_tree(summary.trees[representative_tree_idx],'#' + rep_tree_container_id);
  });
  if(!Config.report_small_clusters){
    var num_unclustered = 0;
    Object.keys(separated_clusters.small).forEach(function(cluster_idx) {
      var C = summary.clusters[cluster_idx];
      num_unclustered = num_unclustered + C.members.length
    })
    var propTrees = (num_unclustered/Object.keys(summary.trees).length).toFixed(2);
    var entries = ["Unclustered"].concat(propTrees).concat("--").map(function(entry) {
      return '<td>' + entry + '</td>';
    }); 
    $('<tr/>').html(entries.join('')).appendTo(cluster_table);
  }
}

ClusterViewer.prototype._draw_rep_tree = function(tree,container_id) {
  var structure = tree.structure;
  var populations = tree.populations;
  var node_ids = Object.keys(populations).map(function(k) {
      return parseInt(k, 10);
  });
  var root_id = Math.min.apply(Math, node_ids);
  var tree_plotter = new TreePlotter();
  var tree_structure = tree_plotter._generate_tree_struct(structure, populations, root_id);
  tree_plotter.draw_tree(tree_structure, container_id, padding = [0, 10, 0, 10], w0 = 150, h0 = 75, radius_scalar = 1/7, include_node_identifier=false);
}