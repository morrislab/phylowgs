function TreePlotter() {
}

TreePlotter.prototype._render_summary_table = function(structure, populations, num_samples, sample_names, root_id) {
  var pop_ids = Util.sort_ints(Object.keys(populations));
  var summary_table = $('#snippets .tree-summary').clone().appendTo('#container');
  summary_table.find('.cellprev').attr('colspan', num_samples);
  summary_table.find('.ccf').attr('colspan', num_samples);
  var header = summary_table.find('thead');
  summary_table = summary_table.find('tbody');

  // Empty cells for Node, SSMs, and CNVs columns.
  var samps_header = ['&mdash;', '&mdash;', '&mdash;'];
  // Do once for cellular prevalence, and once for CCF.
  for(var i = 0; i < 2; i++) {
    sample_names.forEach(function(S) {
      samps_header.push(S);
    });
  }
  var samps_header = samps_header.map(function(entry) {
    return '<th>' + entry + '</th>';
  });
  $('<tr/>').html(samps_header.join('')).appendTo(header);

  var self = this;
  pop_ids.forEach(function(pop_id) {
    var pop = populations[pop_id];
    var cp = pop.cellular_prevalence.map(function(E) { return E.toFixed(3); });
    var ccf = self._calc_ccf(structure, populations, pop_id, root_id);
    ccf = ccf.map(function(E) { return E.toFixed(3); });

    var entries = [pop_id].concat([pop.num_ssms, pop.num_cnvs]).concat(cp).concat(ccf).map(function(entry) {
      return '<td>' + entry + '</td>';
    });
    $('<tr/>').html(entries.join('')).appendTo(summary_table);
  });
}

TreePlotter.prototype._plot_pop_trajectories = function(structure, populations, num_samples, sample_names, hidden_samples, root_id) {
  var pop_ids = Util.sort_ints(Object.keys(populations));
  var num_cancer_pops = pop_ids.length - 1;

  // Don't plot for single-sample data.
  if(num_samples <= 1)
    return;

  var data = new google.visualization.DataTable();
  data.addColumn('string', 'Sample');
  for(var popidx = 1; popidx <= num_cancer_pops; popidx++) {
    data.addColumn('number', 'Population ' + popidx);
  }

  var self = this;
  var ccf = pop_ids.map(function(pop_id) {
    return self._calc_ccf(structure, populations, pop_id, root_id);
  });
  // Remove CCFs for non-cancerous first element, which will always be 1.
  ccf.shift();

  if(num_samples !== sample_names.length) {
    throw 'Improper number of samples provided';
  }
  var data_vals = Util.transpose([sample_names].concat(ccf));
  var visible = [];
  for(var i = 0; i < data_vals.length; i++) {
    var samp = data_vals[i][0];
    if(hidden_samples.indexOf(samp) === -1) {
      visible.push(data_vals[i]);
    } else {
      console.log('Hiding ' + samp);
    }
  }
  data.addRows(visible);

  var container = $('<div/>').appendTo('#container');
  // hAxis.minValue and vAxis.title attributes don't work with Material charts.
  // But the Material charts look so darn sexy, I'm willing to make this
  // sacrifice.
  var options = {
    chart: {
      title: 'Cancer cell fraction trajectories'
    },
    width: container.width(),
    height: 650,
    hAxis: {
      minValue: 1,
      slantedText: true,
      slantedTextAngle: 90,
    },
    vAxis: {
      title: 'Cancer cell fraction',
    },
  };

  //var chart = new google.charts.Line(container.get(0));
  // Uncomment this to switch to using pre-Material charts, which have more
  // options (like, oh, you know, titling the axes) but are less pretty and
  // don't have the hover-over-legend-to-highlight-line function.
  var chart = new google.visualization.LineChart(container.get(0));
  chart.draw(data, options);
}

TreePlotter.prototype.draw = function(populations, structure, root_id, params) {
  if(params !== undefined) { }

  var pop_ids = Util.sort_ints(Object.keys(populations));
  var num_samples = populations[pop_ids[0]].cellular_prevalence.length;

  if(!params || !params.samples) {
    var sample_names = (new Array(num_samples)).fill(0).map(function(val, idx) {
      return 'Sample ' + (idx + 1);
    });
  } else {
    var sample_names = params.samples;
  }
  if(!params || !params.hidden_samples) {
    var hidden_samples = [];
  } else {
    var hidden_samples = params.hidden_samples;
  }

  // This may not exist in the JSON, as it was an extension to the format on
  // June 27, 2017.
  if(typeof root_id === 'undefined') {
    // Find smallest node ID.
    var node_ids = Object.keys(populations).map(function(k) {
        return parseInt(k, 10);
    });
    root_id = Math.min.apply(Math, node_ids);
  }

  var root = this._generate_tree_struct(structure, populations, root_id);
  this._draw_tree(root, '#container');
  this._plot_pop_trajectories(structure, populations, num_samples, sample_names, hidden_samples, root_id);
  this._render_summary_table(structure, populations, num_samples, sample_names, root_id);
}

TreePlotter.prototype._calc_ccf = function(structure, populations, pop_id, root_id) {
  var pop_ids = Util.sort_ints(Object.keys(populations));
  var num_samples = populations[pop_ids[0]].cellular_prevalence.length;

  // CCF of non-cancerous population should be zero.
  if(parseInt(pop_id, 10) === 0)
    return (new Array(num_samples)).fill(0);

  var clonal_pidxs = structure[root_id];
  var purities = [];

  // Don't assume that populations[1] will have the maximum cellular
  // prevalence, as this may not be true for polyclonal tumors.
  for(var sampidx = 0; sampidx < num_samples; sampidx++) {
    var purity = 0;
    clonal_pidxs.forEach(function(pidx) {
        purity += populations[pidx].cellular_prevalence[sampidx];
    });
    purities.push(purity);
  }

  var cps = populations[pop_id].cellular_prevalence;
  var ccf = [];
  for(var sampidx = 0; sampidx < num_samples; sampidx++) {
    ccf.push(cps[sampidx] / purities[sampidx]);
  }

  return ccf;
}

TreePlotter.prototype._draw_tree = function(root, container) {
  // horiz_padding should be set to the maximum radius of a node, so a node
  // drawn on a boundry won't go over the canvas edge. Since max_area = 8000,
  // we have horiz_padding = sqrt(8000 / pi) =~ 51.
  var horiz_padding = 51;
  var m = [10, horiz_padding, 10, horiz_padding],
      w = 800 - m[1] - m[3],
      h = 600 - m[0] - m[2],
      i = 0;

  // Compute the new tree layout.
  var tree = d3.tree().size([h, w]);
  root = tree(d3.hierarchy(root));
  root.descendants().sort(function(a, b) {
    return d3.ascending(parseInt(a.data.name, 10), parseInt(b.data.name, 10));
  });
  var svg = d3.select(container).html('').append('svg:svg')
      .attr('width', w + m[1] + m[3])
      .attr('height', h + m[0] + m[2]);
  var vis = svg.append('svg:g')
      .attr('transform', 'translate(' + m[3] + ',' + m[0] + ')');

  // Update the nodes.
  var node = vis.selectAll('g.node')
      .data(root.descendants(), function(d) { return d.data.name; });

  // Enter any new nodes at the parent's previous position.
  var nodeEnter = node.enter().append('svg:g');
  nodeEnter.attr('class', 'node')
    .attr('transform', function(d) { return 'translate(' + d.y + ',' + d.x + ')'; });
  nodeEnter.append('svg:circle')
      .attr('r', function(d) { return d.data.radius; });
  nodeEnter.append('svg:text')
      .attr('font-size', '30')
      .attr('dominant-baseline', 'central')
      .attr('text-anchor', 'middle')
      .text(function(d) { return d.data.name; });

  // Update the links.
  var link = vis.selectAll('path.link')
      .data(root.links(), function(d) { return d.target.data.name; })
      .attr('stroke-width', '1.5px')

  // Enter any new links at the parent's previous position.
  link.enter().insert('svg:path', 'g')
    .attr('class', 'link')
    .attr('stroke', '#aaa')
    .attr('d', d3.linkHorizontal().x(function(d) {
      return d.y;
    }).y(function(d) {
      return d.x;
    }));
}

TreePlotter.prototype._find_max_ssms = function(populations) {
  var max_ssms = 0;
  for(var pop_id in populations) {
    var pop = populations[pop_id];
    if(pop.num_ssms > max_ssms)
      max_ssms = pop.num_ssms;
  }
  return max_ssms;
}

TreePlotter.prototype._generate_tree_struct = function(adjlist, pops, root_id) {
  var max_ssms = this._find_max_ssms(pops);

  var _add_node = function(node_id, struct) {
    struct.name = node_id;

    var num_ssms = pops[node_id]['num_ssms'];
    struct.radius = TreeUtil.calc_radius(num_ssms /  max_ssms);

    if(typeof adjlist[node_id] === 'undefined') {
      return;
    }
    struct.children = [];
    adjlist[node_id].forEach(function(child_id) {
      var child = {};
      struct.children.push(child);
      _add_node(child_id, child);
    });
  };

  var root = {};
  _add_node(root_id, root);
  return root;
}
