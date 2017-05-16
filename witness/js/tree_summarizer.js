function TreeSummarizer() {
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

TreeSummarizer.prototype._render_lin_idx_vs_branch_idx = function(tree_summary) {
  if(!tree_summary.hasOwnProperty('tree_densities')) {
    return;
  }

  var best_tree_idx = this._find_best_tree(tree_summary.tree_densities);
  var xpoints = [];
  var ypoints = [];
  var marker_symbols = [];
  var marker_sizes = [];
  var labels = [];
  var colours = [];
  var epsilon = 0.000001;

  Object.keys(tree_summary.trees).forEach(function(tidx) {
    tidx = parseInt(tidx, 10);
    var T = tree_summary.trees[tidx];

    labels.push('Tree ' + tidx);
    xpoints.push(T.clustering_index);
    // Epsilon prevents division by zero when CI = 1 (and so BI = LI = 0).
    ypoints.push(T.branching_index / (T.branching_index + T.linearity_index + epsilon));
    marker_symbols.push(tidx === best_tree_idx ? 'cross' : 'dot');
    marker_sizes.push(tidx === best_tree_idx ? 30 : 6);
    colours.push(tree_summary.tree_densities[tidx]);
  });

  var traces = [{
    x: xpoints,
    y: ypoints,
    name: 'points',
    mode: 'markers',
    type: 'scatter',
    text: labels,
    marker: {
      symbol: marker_symbols,
      size: marker_sizes,
      line: { width: 0 },
      colorscale: 'Viridis',
      color: colours,
    },
  }];
  var layout = {
    title: 'Clustering degree vs. branching degree (best tree: ' + best_tree_idx + ')',
    height: 1000,
    xaxis: { title: 'CI'},
    yaxis: { title: 'BI / (BI + LI)'},
    hovermode: 'closest',
    plot_bgcolor: '#440154',
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

    self._render_lin_idx_vs_branch_idx(summary);
    self._render_cell_prevs(cell_prevs);
    self._render_ssm_counts(ssm_counts);
    self._render_pop_counts(pop_counts, min_ssms);
  });
}
