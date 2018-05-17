function TreeSummarizer() {
}

TreeSummarizer.prototype._draw_rep_tree = function(tree,container_id) {
  var structure = tree.structure;
  var populations = tree.populations;
  var node_ids = Object.keys(populations).map(function(k) {
      return parseInt(k, 10);
  });
  var root_id = Math.min.apply(Math, node_ids);
  var tree_plotter = new TreePlotter();
  var tree_structure = tree_plotter._generate_tree_struct(structure, populations, root_id);
  tree_plotter.draw_tree(tree_structure, container_id,padding = [0, 10, 0, 10], w0 = 150, h0 = 75, radius_scalar = 1/7, include_node_identifier=false);
}

TreeSummarizer.prototype._render_vafs = function(dataset) {
  var muts_path = dataset.muts_path;

  d3.json(muts_path, function(muts) {
    var vafs = [];
    for(var ssm_id in muts.ssms) {
      var ssm = muts.ssms[ssm_id];
      var ssm_vafs = [];
      for(var i = 0; i < ssm.ref_reads.length; i++) {
        var a = ssm.ref_reads[i];
        var d = ssm.total_reads[i];
        ssm_vafs.push((d - a)/d);
      }
      vafs.push([Util.mean(ssm_vafs)]);
    }

    var data = new google.visualization.DataTable();
    data.addColumn('number', 'VAF');
    data.addRows(vafs);

    var x_min = 0;
    var x_max = Math.max(1.0, Util.array_max(vafs));
    var container = $('<div/>').prependTo('#container');
    var options = {
      title: 'VAFs (' + vafs.length + ' variants)',
      histogram: { bucketSize: 0.03 },
      fontSize: Config.font_size,
      hAxis: {
        title: 'VAF',
        viewWindow: {
          min: x_min,
          max: x_max
        }
      },
      vAxis: {
        title: 'Number of variants',
      },
      width: container.width(),
      height: 450,
    };

    // Use prependTo() to ensure VAFs are always first plot displayed.
    var chart = new google.visualization.Histogram(container.get(0));
    chart.draw(data, options);
  });
}

TreeSummarizer.prototype._extract_pops_with_top_cell_prevs = function(populations, desired_pops) {
  var pops = [];
  for(var pop_idx in populations) {
    pops.push(populations[pop_idx]);
  }
  pops.sort(function(a, b) {
    // "b - a" means reverse sort.
    return Util.mean(b.cellular_prevalence) - Util.mean(a.cellular_prevalence);
  });

  // Exclude clonal population.
  var sliced = pops.slice(1, desired_pops + 1);
  while(sliced.length < desired_pops)
    sliced.push(null);
  return sliced;
}

TreeSummarizer.prototype._render_cluster_table = function(summary) {
  var cluster_table = $('#cluster-table .tree-summary').clone().appendTo('#container');
  cluster_table = cluster_table.find('tbody');
  var clusters = this._separate_clusters_by_size(summary.clusters, Object.keys(summary.trees).length, Config.tiny_cluster_criteria)
  self = this;
  var clust_num = 0;
  Object.keys(clusters.large).forEach(function(cluster_idx) {
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
  if(Config.group_tiny_clusters & (Object.keys(clusters.small).length != 0)){
    var num_unclustered = 0;
    Object.keys(clusters.small).forEach(function(cluster_idx) {
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

TreeSummarizer.prototype._render_cell_prevs = function(cell_prevs) {
  for(var i = 0; i < cell_prevs.length; i++) {
    var data = new google.visualization.DataTable();
    data.addColumn('number', 'Cellular prevalence');
    data.addRows(cell_prevs[i]);

    var x_min = 0;
    var x_max = 1.0;
    var container = $('<div/>').appendTo('#container');
    var options = {
      title: 'Cellular prevalence (cancerous population ' + (i + 1) + ') (' + cell_prevs[i].length + ' values)',
      fontSize: Config.font_size,
      hAxis: {
        title: 'Cellular prevalence',
      },
      vAxis: {
        title: 'Trees',
      },
      width: container.width(),
      height: 450,
    };

    var chart = new google.visualization.Histogram(container.get(0));
    chart.draw(data, options);
  }
}

TreeSummarizer.prototype._render_ssm_counts = function(ssm_counts) {
  for(var i = 0; i < ssm_counts.length; i++) {
    var data = new google.visualization.DataTable();
    data.addColumn('number', 'SSMs');
    data.addRows(ssm_counts[i]);

    var container = $('<div/>').appendTo('#container');
    var options = {
      title: 'Number of SSMs (cancerous population ' + (i + 1) + ') (' + ssm_counts[i].length + ' values)',
      fontSize: Config.font_size,
      hAxis: {
        title: 'SSMs',
      },
      vAxis: {
        title: 'Trees',
      },
      width: container.width(),
      height: 450,
    };

    var chart = new google.visualization.Histogram(container.get(0));
    chart.draw(data, options);
  }
}

TreeSummarizer.prototype._render_pop_counts = function(pop_counts, min_ssms) {
  var histogram = {};
  var min_count = pop_counts.length, max_count = 0;
  pop_counts.forEach(function(count) {
    if(count < min_count)
      min_count = count;
    if(count > max_count)
      max_count = count;
    if(count in histogram) {
      histogram[count]++;
    } else {
      histogram[count] = 1;
    }
  });

  var rows = [];
  for(var i = min_count; i <= max_count; i++) {
    if(i in histogram)
      rows.push([i.toString(), histogram[i]]);
    else
      rows.push([i.toString(), 0]);
  }

  var data = new google.visualization.DataTable();
  data.addColumn('string', 'Populations');
  data.addColumn('number', 'Count');
  data.addRows(rows);

  var container = $('<div/>').appendTo('#container');
  var options = {
    title: 'Distribution of cancerous populations (' + pop_counts.length + ' values)',
    fontSize: Config.font_size,
    hAxis: {
      title: 'Number of cancerous populations',
    },
    vAxis: {
      title: 'Trees',
    },
    width: container.width(),
    height: 450,
  };

  var chart = new google.visualization.ColumnChart(container.get(0));
  chart.draw(data, options);
}

TreeSummarizer.prototype._extract_indices = function(tree_summ) {
  this.linearity_indices = {};
  this.branching_indices = {};
  this.clustering_indices = {};

  var self = this;
  Object.keys(tree_summ).forEach(function(tidx) {
    self.linearity_indices[tidx] = tree_summ[tidx].linearity_index;
    self.branching_indices[tidx] = tree_summ[tidx].branching_index;
    self.clustering_indices[tidx] = tree_summ[tidx].clustering_index;
  });
  this.num_trees = Object.keys(this.linearity_indices).length;
}

TreeSummarizer.prototype._determine_cluster_colour_ordering = function(clusters,colors,xData,yData){
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

TreeSummarizer.prototype._separate_clusters_by_size = function(clusters, total_num_trees, cluster_size_criteria){
  if(!Config.group_tiny_clusters){
    return {small: [], large: clusters};
  }
  var small_clusters = {};
  var large_clusters = {};
  Object.keys(clusters).forEach(function(cidx){
    var C = clusters[cidx];
    if(C.members.length < total_num_trees*cluster_size_criteria){
      small_clusters[cidx] = C;
    }else{
      large_clusters[cidx] = C;
    }
  })
  return {small: small_clusters, large: large_clusters};
}

TreeSummarizer.prototype._find_best_tree = function(densities) {
  var max_density = 0;
  var best_tidx = null;
  Object.keys(densities).forEach(function(tidx) {
    var density = densities[tidx];
    if(density > max_density) {
      max_density = density;
      best_tidx = tidx;
    }
  });
  if(best_tidx == null) {
    throw "best_tidx is null";
  }
  return parseInt(best_tidx, 10);
}

TreeSummarizer.prototype._calc_xy_data = function(trees){
  var xData = [];
  var yData = [];
  var epsilon = 0.000001
  Object.keys(trees).forEach(function(tidx){
    var T = trees[tidx];
    xData.push(T.clustering_index)
    // Epsilon prevents division by zero when CI = 1 (and so BI = LI = 0)
    yData.push(T.branching_index / (T.branching_index + T.linearity_index + epsilon));
  });
  return {x:xData,y:yData}
}

TreeSummarizer.prototype._create_tiny_cluster_scatter_trace = function(clusters, xData, yData) {
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

TreeSummarizer.prototype._create_cluster_scatter_traces = function(clusters, cluster_colours, xData, yData) {
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
      marker_symbols.push(rep_tree_index === tidx ? 'cross' : 'dot');
      marker_sizes.push(rep_tree_index === tidx ? 30 : 6);
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
      showlegend: !Config.show_tree_cluster_contours,
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

TreeSummarizer.prototype._create_cluster_contour_traces = function(clusters, cluster_colours) {

  var traces = [];
  var clust_count = 0;
  var colour_idx = 0;
  Object.keys(clusters).forEach(function(cidx) {
    //cidx = parseInt(cidx, 10);
    clust_count = clust_count + 1;
    var C = clusters[cidx];
    var mean = C.ellipse.mean //Center of the ellipse
    var angle = C.ellipse.angle //angle of the ellipse wrt the x axis.
    var maj_axis = C.ellipse.major_axis //the distance between the center of the ellipse and the farthest point, ie, the highest variance of the gaussian
    var min_axis = C.ellipse.minor_axis //the distance between the center of the ellipse and the closest point, ie, the lowest variance of the gaussian
    var xpoints = [];
    var ypoints = [];
    for(var theta = 0.; theta <= 2.*Math.PI + 0.05; theta+=2.*Math.PI/360.){
      //equations for x and y values of the ellipse described by the mean, angle, major axis and minor axis, as described as a function of theta.
      xpoints.push(mean[0] + maj_axis*Math.cos(angle)*Math.cos(theta) - min_axis*Math.sin(angle)*Math.sin(theta))
      ypoints.push(mean[1] + maj_axis*Math.sin(angle)*Math.cos(theta) + min_axis*Math.cos(angle)*Math.sin(theta))
    }
    //Descriptor of traces for now. May update depending on visualization requirements.
    traces.push({
      x: xpoints,
      y: ypoints,
      name: 'cluster ' + clust_count,
      text: 'cluster ' + clust_count,
      type: 'scatter',
      mode: 'lines',
      visible: Config.show_tree_cluster_contours,
      line: {color: cluster_colours[colour_idx]}
    })
    colour_idx = colour_idx + 1;
  });
  return traces
}

TreeSummarizer.prototype._is_using_old_data = function(summary){
    // Check to see if using older data. If so, may not be able to plot some graphs.
    var using_old_data = false;
    if(!summary.hasOwnProperty('clusters')){
      using_old_data = true;
    };
    Object.keys(summary.clusters).forEach(function(cKey){
      if (!summary.clusters[cKey].hasOwnProperty('members')){
        using_old_data = true;
        return
      }
    })
    return using_old_data;
}

TreeSummarizer.prototype._render_lin_idx_vs_branch_idx = function(tree_summary) {
  if(!tree_summary.hasOwnProperty('tree_densities') | !tree_summary.hasOwnProperty('clusters')) {
    return;
  }

  // Calculcate the xData and yData from the indicies
  var xyData = this._calc_xy_data(tree_summary.trees);

  var best_tree_idx = this._find_best_tree(tree_summary.tree_densities);
  var clusters = this._separate_clusters_by_size(tree_summary.clusters, Object.keys(tree_summary.trees).length, Config.tiny_cluster_criteria);
  var cluster_colours = this._determine_cluster_colour_ordering(clusters.large, Config.cluster_colours, xyData.x, xyData.y);
  // Determine the ellipse traces that define the cluster contours
  var ellipse_traces = this._create_cluster_contour_traces(clusters.large, cluster_colours);
  // Create the traces for the data points for clustered trees
  var scatter_traces = this._create_cluster_scatter_traces(clusters.large, cluster_colours, xyData.x, xyData.y);
  // Create the trace for the data points for unclustered trees
  scatter_traces.push(this._create_tiny_cluster_scatter_trace(clusters.small, xyData.x, xyData.y))
  var traces = ellipse_traces.concat(scatter_traces);
  
  //Use the traces and plotting options to actually plot our data.
  var layout = {
    title: "Clustering Degree vs. Branching Degree (best tree: " + best_tree_idx + ")",
    height: 1000,
    xaxis: { title: 'CI'},
    yaxis: { title: 'BI/(LI+BI)'},
    hovermode: 'closest',
    plot_bgcolor: '#440154',
    showlegend: true
  };
  var container = document.querySelector('#container');
  var plot_container = document.createElement('div');
  container.appendChild(plot_container);
  Plotly.newPlot(plot_container, traces, layout);
}

TreeSummarizer.prototype.render = function(dataset) {
  this._render_vafs(dataset);

  var pops_to_examine = 3;
  var min_ssms = 3;

  var pop_counts = [];
  var cell_prevs = new Array(pops_to_examine);
  var ssm_counts = new Array(pops_to_examine);
  var lin_idx = [];
  var branch_idx = [];
  var tree_idx = [];

  for(var i = 0; i < pops_to_examine; i++) {
    cell_prevs[i] = [];
    ssm_counts[i] = [];
  }

  var self = this;
  d3.json(dataset.summary_path, function(summary) {
    Object.keys(summary.trees).forEach(function(tidx) {
      var populations = summary.trees[tidx].populations;
      var num_pops = 0;
      for(var pop_idx in populations) {
        var pop = populations[pop_idx];
        if(pop.num_cnvs > 0 || pop.num_ssms >= min_ssms) {
          num_pops++;
        }
      }
      pop_counts.push(num_pops);

      var pops = self._extract_pops_with_top_cell_prevs(populations, pops_to_examine);
      for(var i = 0; i < pops.length; i++) {
        if(pops[i] !== null) {
          cell_prevs[i].push([Util.mean(pops[i].cellular_prevalence)]);
          ssm_counts[i].push([pops[i].num_ssms]);
        }
      }
    });
    if(!self._is_using_old_data(summary)){
      // These functions require that the summ file be created with the latest version of write_results.
      // If not, then just skip these sections.
      self._render_cluster_table(summary);
      self._render_lin_idx_vs_branch_idx(summary);
    }
    self._render_cell_prevs(cell_prevs);
    self._render_ssm_counts(ssm_counts);
    self._render_pop_counts(pop_counts, min_ssms);
  });
}
