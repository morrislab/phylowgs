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

TreeSummarizer.prototype._cov_matrix = function(X) {
  var XT = transposeMatrix(X);
  var N = X.length;
  var D = XT.length;

  var means = [];
  for(var i = 0; i < D; i++) {
    var sum = XT[i].reduce(function(a, b) { return a + b; });
    means.push(sum / N);
  }
  var residuals = [];
  for(var i = 0; i < D; i++) {
    residuals[i] = [];
    for(var j = 0; j < N; j++) {
      residuals[i].push(XT[i][j] - means[i]);
    }
  }

  var cov = zeros(D, D);
  for(var i = 0; i < D; i++) {
    for(var j = 0; j < D; j++) {
      var res_sum = 0;
      for(var k = 0; k < N; k++) {
        res_sum += (XT[i][k] - means[i]) * (XT[j][k] - means[j]);
      }
      cov[i][j] = res_sum / (N - 1);
    }
  }

  return cov;
}

TreeSummarizer.prototype._gaussian_pdf = function(X, mu, sigma) {
  var lambda = inverseMatrix(sigma);
  var detsigma = determinant(sigma);

  var density = X.map(function(x) {
    var D = x.length;
    var R = [x.map(function(x, i) { return x - mu[i]; })]; // 1xD
    var dens = Math.exp(-0.5 * matrixMultiply(R, matrixMultiply(lambda, transposeMatrix(R))));
    dens /= Math.sqrt(Math.pow(2 * Math.PI, D) * detsigma);
    return dens;
  });
  return density;
}

TreeSummarizer.prototype._zerovec = function(N) {
  return Array.apply(null, Array(N)).map(Number.prototype.valueOf,0);
}

TreeSummarizer.prototype._kde = function(points) {
  var N = points.length;
  var D = points[0].length;

  // Scott's Rule from SciPy
  var bandwidth = Math.pow(N, -1/(D + 4));
  var sigma = this._cov_matrix(points);
  sigma = matrixScalarMultiply(sigma, Math.pow(bandwidth, 2));

  var density = this._zerovec(N);
  var self = this;
  points.forEach(function(P) {
    var density_P = self._gaussian_pdf(points, P, sigma);
    for(var i = 0; i < density_P.length; i++) {
      density[i] += density_P[i];
    }
  });
  return density;
}

TreeSummarizer.prototype._render_lin_idx_vs_branch_idx = function(lin_idx, branch_idx, tree_idx) {
  var marker_symbols = [];
  var marker_sizes = [];

  var btf = new BestTreeFinder(lin_idx, branch_idx, tree_idx);
  var indices = btf.make_indices();
  var mean_index = btf.calc_mean_index();
  var best_tree_idx = btf.find_best_tree();

  var labels = tree_idx.map(function(T) { return 'Tree ' + T; });

  for(var i = 0; i < lin_idx.length; i++) {
    marker_symbols.push(i === best_tree_idx ? 'cross' : 'dot');
    marker_sizes.push(i === best_tree_idx ? 30 : 6);
  }

  var xpoints = lin_idx.slice(); // Duplicate array
  var ypoints = branch_idx.slice();
  xpoints.push(mean_index[0]);
  ypoints.push(mean_index[1]);
  marker_symbols.push('diamond');
  marker_sizes.push(30);
  labels.push('Mean');

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
      color: '#ff5858',
    },
  },
  {
    x: xpoints,
    y: ypoints,
    ncontours: 20,
    colorscale: 'Viridis',
    type: 'histogram2dcontour',
  }];
  var layout = {
    title: 'Branching index vs. linearity index (best tree: ' + best_tree_idx + ')',
    height: 1000,
    xaxis: { title: 'Linearity index'},
    yaxis: { title: 'Branching index'},
    hovermode: 'closest',
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
    var tidxs = Object.keys(summary.trees);
    tidxs.sort(function(a, b) { return parseInt(a, 10) - parseInt(b, 10); });

    tidxs.forEach(function(tidx) {
      lin_idx.push(summary.trees[tidx].linearity_index);
      branch_idx.push(summary.trees[tidx].branching_index);
      tree_idx.push(tidx);

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

    self._render_lin_idx_vs_branch_idx(lin_idx, branch_idx, tree_idx);
    self._render_cell_prevs(cell_prevs);
    self._render_ssm_counts(ssm_counts);
    self._render_pop_counts(pop_counts, min_ssms);
  });
}
