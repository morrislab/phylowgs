Config = {
  font_size: 15
};

var MARKER = '<defs><marker id="head" orient="auto" markerWidth="2" markerHeight="4" refX="0" refY="2"><path d="M0,0 V4 L2,2 Z" fill="#000" /></marker></defs>';

function mean(list) {
  return list.reduce(function(a, b) { return a + b; }) / list.length;
}

function TreePlotter() {
}

TreePlotter.prototype._sort_numeric = function(arr) {
  return arr.sort(function(a, b) {
    return parseInt(a, 10) - parseInt(b, 10);
  });
}

TreePlotter.prototype._calc_ccf = function(tree, pop_id) {
  if(parseInt(pop_id, 10) === 0)
    return 0;

  // This only works for monoclonal trees, but it's the desired behaviour according to Quaid.
  var cellularity = mean(tree.populations[1].cellular_prevalence);
  var ccf = mean(tree.populations[pop_id].cellular_prevalence) / cellularity;
  return ccf;
}

TreePlotter.prototype.render = function(dataset) {
  $('#tree-list').show();
  var tree_container = $('#trees tbody');

  var tplotter = this;
  d3.json(dataset.summary_path, function(summary) {
    var tree_indices = tplotter._sort_numeric(Object.keys(summary.trees));
    tree_container.empty();
    tree_indices.forEach(function(tidx) {
      var row = '<td class="tree-index">' + tidx + '</td>'
        + '<td class="tree-llh">' + summary.trees[tidx].llh.toFixed(1) + '</td>'
        + '<td class="tree-nodes">' + Object.keys(summary.trees[tidx].populations).length + '</td>';
      $('<tr/>').html(row).appendTo(tree_container);
    });

    $('#trees').stupidtable();

    var already_autosorted = false;
    $('#trees').bind('aftertablesort', function() {
      if(already_autosorted)
        return;
      tree_container.find('tr:first').click();
      already_autosorted = true;
    });

    // If direction not specified, this can end up being ascending or
    // descending sort, depending on prior sort state of table.
    $('#tree-llh').stupidsort('desc');

    tree_container.find('tr').click(function(evt) {
      evt.preventDefault();
      var self = $(this);
      self.siblings().removeClass('active');
      self.addClass('active');

      var tidx = self.find('.tree-index').text();
      var root = tplotter._generate_tree_struct(summary.trees[tidx].structure, summary.trees[tidx].populations);
      tplotter._draw_tree(root);

      var summary_table = $('#snippets .tree-summary').clone().appendTo('#container').find('tbody');
      var pop_ids = tplotter._sort_numeric(Object.keys(summary.trees[tidx].populations));
      pop_ids.forEach(function(pop_id) {
        var pop = summary.trees[tidx].populations[pop_id];
        var ccf = tplotter._calc_ccf(summary.trees[tidx], pop_id);
        var entries = [pop_id, mean(pop.cellular_prevalence).toFixed(3), ccf.toFixed(3), pop.num_ssms, pop.num_cnvs].map(function(entry) {
          return '<td>' + entry + '</td>';
        });
        $('<tr/>').html(entries.join('')).appendTo(summary_table);
      });
    });
    $('#tree-list').scrollTop(0);
  });
}

TreePlotter.prototype._draw_tree = function(root) {
  // horiz_padding should be set to the maximum radius of a node, so a node
  // drawn on a boundry won't go over the canvas edge. Since max_area = 8000,
  // we have horiz_padding = sqrt(8000 / pi) =~ 51.
  var horiz_padding = 51;
  var m = [10, horiz_padding, 10, horiz_padding],
      w = 800 - m[1] - m[3],
      h = 600 - m[0] - m[2],
      i = 0;

  var tree = d3.layout.tree()
      .size([h, w])
      .sort(function(a, b) {
        return d3.ascending(parseInt(a.name, 10), parseInt(b.name, 10));
      });

  var diagonal = d3.svg.diagonal()
      .projection(function(d) { return [d.y, d.x]; });

  var vis = d3.select('#container').html('').append('svg:svg')
      .attr('width', w + m[1] + m[3])
      .attr('height', h + m[0] + m[2])
      .append('svg:g')
      .attr('transform', 'translate(' + m[3] + ',' + m[0] + ')');

  // Compute the new tree layout.
  var nodes = tree.nodes(root);

  // Update the nodes…
  var node = vis.selectAll('g.node')
      .data(nodes, function(d) { return d.name; });

  // Enter any new nodes at the parent's previous position.
  var nodeEnter = node.enter().append('svg:g')
      .attr('class', 'node');

  nodeEnter.append('svg:circle')
      .attr('r', function(d) { return d.radius; });

  nodeEnter.append('svg:text')
      .attr('font-size', '30')
      .attr('dominant-baseline', 'central')
      .attr('text-anchor', 'middle')
      .text(function(d) { return d.name; });

  // Transition nodes to their new position.
  var nodeUpdate = node.attr('transform', function(d) { return 'translate(' + d.y + ',' + d.x + ')'; });
  var nodeExit = node.exit().remove();

  // Update the links…
  var link = vis.selectAll('path.link')
      .data(tree.links(nodes), function(d) { return d.target.name; })
      .attr('stroke-width', '1.5px');

  // Enter any new links at the parent's previous position.
  link.enter().insert('svg:path', 'g')
      .attr('class', 'link')
      .attr('stroke', '#aaa');

  // Transition links to their new position.
  link.attr('d', diagonal);

  // Transition exiting nodes to the parent's new position.
  link.exit().remove();
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

TreePlotter.prototype._generate_tree_struct = function(adjlist, pops) {
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
  var node_ids = Object.keys(pops).map(function(k) {
      return parseInt(k, 10);
  });
  // Find smallest node ID.
  var initial_node = Math.min.apply(Math, node_ids);

  _add_node(initial_node, root);
  return root;
}

function TreeUtil() {
}

TreeUtil.calc_radius = function(scale) {
  var min_area = 700, max_area = 8000;
  var area = Util.calc_in_range(min_area, max_area, scale);
  return Math.sqrt(area / Math.PI);
}

function ClusterPlotter() {
  // horiz_padding should be set to the maximum radius of a node, so a node
  // drawn on a boundry won't go over the canvas edge. Since max_area = 8000,
  // we have horiz_padding = sqrt(8000 / pi) =~ 51.
  var horiz_padding = 51;
  this._M = [10, horiz_padding, 10, horiz_padding],
      this._width = 800 - this._M[1] - this._M[3],
      this._height = 600 - this._M[0] - this._M[2];
}

ClusterPlotter.prototype.render = function(dataset) {
  if(typeof dataset.clusters_path === 'undefined')
    return;

  var cplotter = this;
  d3.json(dataset.clusters_path, function(formatted_clusters) {
    var cluster_indices = Object.keys(formatted_clusters);

    $('#cluster-list').show();
    var cluster_table = $('#clusters tbody');
    cluster_table.empty();

    cluster_indices.forEach(function(cidx) {
      var row = '<td class="cluster-index">' + cidx + '</td>'
        + '<td class="cluster-size">' + formatted_clusters[cidx].members.length + '</td>';
      $('<tr/>').html(row).appendTo(cluster_table);
    });

    $('#clusters').stupidtable();

    var already_autosorted = false;
    $('#clusters').bind('aftertablesort', function() {
      if(already_autosorted)
        return;
      cluster_table.find('tr:first').click();
      already_autosorted = true;
    });

    // If direction not specified, this can end up being ascending or
    // descending sort, depending on prior sort state of table.
    $('#cluster-idx').stupidsort('asc');

    cluster_table.find('tr').click(function(evt) {
      evt.preventDefault();
      var self = $(this);
      self.siblings().removeClass('active');
      self.addClass('active');

      var cidx = self.find('.cluster-index').text();
      cplotter._show_cluster(formatted_clusters[cidx]);
    });
    $('#cluster-list').scrollTop(0);
  });
}

ClusterPlotter.prototype._empty_plot_container = function() {
  if($('#container #plots').length === 0) {
    $('<div/>').attr('id', 'plots').appendTo('#container');
  } else {
    $('#plots').empty();
  }
}

ClusterPlotter.prototype._draw_cluster_plots = function(cluster) {
  this._empty_plot_container();

  var simils = [];
  cluster.members.forEach(function(member) {
    var simil = member[1];
    simils.push(simil);
  });

  this._draw_histogram(simils, 'Tree similarities', 'Similarity', 0, 1);
}

ClusterPlotter.prototype._draw_pop_plots = function(pop) {
  // Don't draw plots for root node, as its parameter values are always fixed.
  if(pop.name === 0)
    return;

  this._empty_plot_container();

  if(!pop.is_bastard) {
    this._draw_histogram(pop.jaccard_scores, 'Jaccard scores (population ' + pop.name + ')', 'Jaccard scores', 0, 1);
  }
  this._draw_histogram(pop.num_ssms, 'Number of SSMs (population ' + pop.name + ')', 'Number of SSMs');

  var mean_cell_prevs = [];
  pop.cellular_prevalences.forEach(function(cp) {
    mean_cell_prevs.push(Util.mean(cp));
  });
  this._draw_histogram(mean_cell_prevs, 'Cellular prevalences (population ' + pop.name + ')', 'Cellular prevalence', 0, 1);
}

ClusterPlotter.prototype._draw_histogram = function(vals, title, xlabel, xmin, xmax) {
  var formatted = [];
  vals.forEach(function(val) {
    formatted.push([val]);
  });

  var data = new google.visualization.DataTable();
  data.addColumn('number', xlabel);
  data.addRows(formatted);

  var view_window = {};
  if(typeof xmin !== 'undefined')
    view_window.min = xmin;
  if(typeof xmax !== 'undefined')
    view_window.max = xmax;

  var options = {
    title: title + ' (' + vals.length + ' values)',
    fontSize: Config.font_size,
    hAxis: {
      title: title,
      viewWindow: view_window,
    },
    vAxis: {
      title: 'Trees',
    },
    width: '100%',
    height: 450,
  };

  var chart = new google.visualization.Histogram($('<div/>').appendTo('#plots').get(0));
  chart.draw(data, options);
}

ClusterPlotter.prototype._show_cluster = function(cluster) {
  var trees_in_cluster = cluster.members.length;

  var max_ssms = 0;
  for(var popidx in cluster.populations) {
    var mean_ssms = Util.mean(cluster.populations[popidx].num_ssms);
    if(mean_ssms > max_ssms)
      max_ssms = mean_ssms;
  }

  var pops = [];
  var self = this;
  var Y = 40;
  Object.keys(cluster.populations).sort().forEach(function(popidx) {
    popidx = parseInt(popidx, 10);
    if(popidx !== pops.length) {
      throw "Expected " + popidx + " populations, but have " + pops.length;
    }

    var popdata = cluster.populations[popidx];
    var trees_with_pop = popdata.trees_with_pop;
    var mean_ssms = Util.mean(popdata.num_ssms);

    var radius = TreeUtil.calc_radius(mean_ssms / max_ssms);
    var pop = {
      x: self._width / 2 - radius,
      y: Y,
      name: popidx,
      radius: radius,
      opacity: trees_with_pop / trees_in_cluster,
      is_bastard: false,
      jaccard_scores: popdata.jaccard_scores,
      num_ssms: popdata.num_ssms,
      cellular_prevalences: popdata.cellular_prevalences,
      children: popdata.children || [],
      bastards: popdata.bastard_subtrees || null
    };
    if(popidx == 0) {
      pop.fixed = true;
    }
    Y += 200;
    pops.push(pop);
  });

  var links = [];
  pops.forEach(function(pop) {
    // Add edges from pops to pops.
    Object.keys(pop.children).forEach(function(childidx) {
      childidx = parseInt(childidx, 10);
      var child = pops[childidx];
      var trees_with_edge = pop.children[childidx];

      links.push({
        source: pop,
        target: child,
        width: Util.calc_in_range(2, 15, trees_with_edge / trees_in_cluster),
        name: pop.name + '_' + child.name
      });
    });

    // Add bastards.
    if(pop.bastards) {
      pop.bastards.forEach(function(bastard) {
        var mean_ssms = Util.mean(bastard.num_ssms);
        var trees_with_bastard = bastard.num_ssms.length;

        var bastard_node = {
          // Assign numerical index following all populations.
          name: pops.length,
          radius: TreeUtil.calc_radius(mean_ssms / max_ssms),
          opacity: trees_with_bastard / trees_in_cluster,
          num_ssms: bastard.num_ssms,
          cellular_prevalences: bastard.cellular_prevalences,
          is_bastard: true
        };
        pops.push(bastard_node);

        links.push({
          source: pop,
          target: bastard_node,
          width: Util.calc_in_range(2, 15, trees_with_bastard / trees_in_cluster),
          name: pop.name + '_' + bastard_node.name
        });
      });
    }
  });

  this._draw(pops, links);
  this._draw_cluster_plots(cluster);
}

ClusterPlotter.prototype._tick = function(link, node, alpha, redraw) {
  link.each(function(d) {
    // Given each directed edge, push the parent up and the child down,
    // resulting in a more tree-like display.
    var K = 10 * alpha;
    // Effect disabled for the moment because it results in nodes overlapping
    // their arrowheads.
    //d.source.y -= K;
    //d.target.y += K;
  });

  var _shorten_for_arrowhead = function(d) {
    var vec = [ d.target.x - d.source.x, d.target.y - d.source.y ];
    var len = Math.sqrt(vec[0]*vec[0] + vec[1]*vec[1]);
    // In early simulation stages, source and target may be at the same
    // position, yielding 0 length.
    if(len === 0) {
      return { x: 0, y: 0 };
    }

    var radius = d.target.radius;
    // Take 2*d.width since markerWidth = 2.
    var offset = radius + 2*d.width;
    // Don't allow lengths < 2*radius.
    var new_len = Math.max(2*radius, len - offset);

    vec[0] *= (new_len / len);
    vec[1] *= (new_len / len);

    return {
      x: d.source.x + vec[0],
      y: d.source.y + vec[1]
    };
  };

  node.attr('transform', function(d) { return 'translate(' + d.x + ',' + d.y + ')'; })
  link.attr('x1', function(d) { return d.source.x; })
    .attr('y1', function(d) { return d.source.y; })
    .attr('x2', function(d) { return _shorten_for_arrowhead(d).x; })
    .attr('y2', function(d) { return _shorten_for_arrowhead(d).y; });
}

ClusterPlotter.prototype._draw = function(pops, links) {
  var self = this;
  var force = d3.layout.force()
      .size([this._width, this._height])
      .charge(-3500)
      .gravity(0.2)
      .linkStrength(0.01);

  var vis = d3.select('#container').html('').append('svg:svg')
      .html(MARKER)
      .attr('width', this._width + this._M[1] + this._M[3])
      .attr('height', this._height + this._M[0] + this._M[2])
      .append('svg:g')
      .attr('transform', 'translate(' + this._M[3] + ',' + this._M[0] + ')');

  force.nodes(pops)
    .links(links)
    .start();

  var self = this;
  // Update the nodes
  var node = vis.selectAll('.node')
      .data(pops, function(d) { return d.name; });
  // Enter any new nodes at the parent's previous position.
  var nodeEnter = node.enter().append('svg:g')
      .classed('node', true)
      .classed('bastard-subtree', function(d) { return d.is_bastard; })
      .classed('population', function(d) { return !d.is_bastard; })
      .attr('opacity', function(d) { return d.opacity; })
      .call(force.drag)
      .on('click', function(pop) {
        // Suppress event if it resulted from a drag-drop event, which is used
        // to reposition the nodes on the canvas.
        if (d3.event.defaultPrevented) return;
        self._draw_pop_plots(pop);
      });
  nodeEnter.append('svg:circle')
      .attr('r', function(d) { return d.radius; });
  nodeEnter.append('svg:text')
      .attr('font-size', '30')
      .attr('dominant-baseline', 'central')
      .attr('text-anchor', 'middle')
      .text(function(d) { return d.name; });
  var nodeExit = node.exit().remove();

  // Update the links
  var link = vis.selectAll('.link')
      .data(links, function(d) { return d.name; });
  link.exit().remove();

  link.enter().insert('svg:line', '.node')
      .attr('class', 'link')
      .attr('stroke', '#000')
      .attr('opacity', 0.3)
      .attr('stroke-width', function(d) { return d.width + 'px'; })
      .attr('marker-end', 'url(#head)');

  // Advance the simulation to achieve a near-stable solution before rendering graph.
  for(var i = 0; i < 200 && force.alpha() > 0.03; i++) {
    force.tick();
  }
  force.on('tick', function(evt) { self._tick(link, node, evt.alpha); });
}

function TreeSummarizer() {
}

function get_url_param(name) {
  name = name.replace(/[\[]/, "\\[").replace(/[\]]/, "\\]");
  var regex = new RegExp("[\\?&]" + name + "=([^&#]*)"),
      results = regex.exec(location.search);
  return results === null ? "" : decodeURIComponent(results[1].replace(/\+/g, " "));
}

function make_parent_active(elem) {
  var parent = elem.parent();
  parent.siblings('li').removeClass('active');
  parent.addClass('active');
}

function Util() {
}

Util.mean = function(arr) {
  var sum = 0;
  for(var i = 0; i < arr.length; i++)
    sum += arr[i];
  return sum / arr.length;
}

Util.array_max = function(arr) {
  return Math.max.apply(null, arr);
}

// scale must be in [0, 1].
Util.calc_in_range = function(min, max, scale) {
  return min + scale*(max - min);
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
      width: 1000,
      height: 450,
    };

    // Use prependTo() to ensure VAFs are always first plot displayed.
    var container = $('<div/>').prependTo('#container').get(0);
    var chart = new google.visualization.Histogram(container);
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
    return mean(b.cellular_prevalence) - mean(a.cellular_prevalence);
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
    var options = {
      title: 'Cellular prevalence (subclone ' + (i + 1) + ') (' + cell_prevs[i].length + ' values)',
      fontSize: Config.font_size,
      hAxis: {
        title: 'Cellular prevalence',
      },
      vAxis: {
        title: 'Trees',
      },
      width: 1000,
      height: 450,
    };

    var container = $('<div/>').appendTo('#container').get(0);
    var chart = new google.visualization.Histogram(container);
    chart.draw(data, options);
  }
}

TreeSummarizer.prototype._render_ssm_counts = function(ssm_counts) {
  for(var i = 0; i < ssm_counts.length; i++) {
    var data = new google.visualization.DataTable();
    data.addColumn('number', 'SSMs');
    data.addRows(ssm_counts[i]);

    var options = {
      title: 'Number of SSMs (subclone ' + (i + 1) + ') (' + ssm_counts[i].length + ' values)',
      fontSize: Config.font_size,
      hAxis: {
        title: 'SSMs',
      },
      vAxis: {
        title: 'Trees',
      },
      width: 1000,
      height: 450,
    };

    var container = $('<div/>').appendTo('#container').get(0);
    var chart = new google.visualization.Histogram(container);
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

  var options = {
    title: 'Distribution of cancerous populations (' + pop_counts.length + ' values)',
    fontSize: Config.font_size,
    hAxis: {
      title: 'Number of cancerous populations',
    },
    vAxis: {
      title: 'Trees',
    },
    width: 1000,
    height: 450,
  };

  var container = $('<div/>').appendTo('#container').get(0);
  var chart = new google.visualization.ColumnChart(container);
  chart.draw(data, options);
}

TreeSummarizer.prototype.render = function(dataset) {
  this._render_vafs(dataset);

  var pops_to_examine = 3;
  var min_ssms = 3;

  var pop_counts = [];
  var cell_prevs = new Array(pops_to_examine);
  var ssm_counts = new Array(pops_to_examine);
  for(var i = 0; i < pops_to_examine; i++) {
    cell_prevs[i] = [];
    ssm_counts[i] = [];
  }

  var self = this;
  d3.json(dataset.summary_path, function(summary) {
    for(var tidx in summary.trees) {
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
          cell_prevs[i].push([mean(pops[i].cellular_prevalence)]);
          ssm_counts[i].push([pops[i].num_ssms]);
        }
      }
    }

    self._render_cell_prevs(cell_prevs);
    self._render_ssm_counts(ssm_counts);
    self._render_pop_counts(pop_counts, min_ssms);
  });
}

function Interface() {
  this._activate_filters();
  this._activate_navbar();
  this._load_samples();

  this._available_renderers = {
    'summarizer': new TreeSummarizer(),
    'tree_plotter': new TreePlotter(),
    'cluster_plotter': new ClusterPlotter()
  };
}

Interface.prototype._activate_filters = function() {
  $('.filter').keyup(function(evt) {
    var self = $(this);
    var filter_text = self.val();
    var elems = self.parents('.sidebar').find('.nav-sidebar').find('li');
    elems.hide();
    elems.filter(function() {
      return !($(this).text().indexOf(filter_text) === -1);
    }).show();
  });
}

Interface.prototype._activate_navbar = function() {
  var iface = this;
  $('.navbar-nav a').click(function(evt) {
    evt.preventDefault();
    var self = $(this);
    make_parent_active(self);

    var mapping = {
      'nav-tree-summaries': 'summarizer',
      'nav-tree-viewer': 'tree_plotter',
      'nav-clustered-trees': 'cluster_plotter'
    };
    iface._renderer = null;
    for(var nav_class in mapping) {
      if(self.hasClass(nav_class)) {
        iface._renderer = mapping[nav_class];
        break;
      }
    }

    iface._render();
  });
}

Interface.prototype._render = function() {
  if(this._renderer !== null && this._dataset !== null) {
    $('#container').empty();
    $('.secondary-sidebar').hide();
    this._available_renderers[this._renderer].render(this._dataset);
  }
}

Interface.prototype._load_samples = function() {
  var run_container = $('#runs');
  var sample_container = $('#samples');
  //
  // Show tree summaries by default.
  this._renderer = 'summarizer';
  this._dataset = null;

  var iface = this;
  d3.json('data/index.json', function(data_index) {
    Object.keys(data_index).sort().forEach(function(run_name) {
      var li = $('<li/>').appendTo(run_container);
      $('<a/>').text(run_name).attr('href', '#').appendTo(li)
    });

    run_container.find('a').click(function(evt) {
      evt.preventDefault();
      var self = $(this);
      make_parent_active(self);
      var run_name = self.text();

      var prev_dataset = iface._dataset;
      iface._dataset = null;
      sample_container.empty();
      // Remove any visible results from another run.
      $('#container').empty();
      $('.page-header').text(run_name);

      data_index[run_name].sort(function(a, b) {
        if(a.name > b.name) return 1;
        if(a.name < b.name) return -1;
        return 0;
      }).forEach(function(ds) {
        var li = $('<li/>').appendTo(sample_container);
        $('<a/>').text(ds.name).attr('href', '#').data('dataset', ds).appendTo(li)
      });

      sample_container.find('a').click(function(evt) {
        var self = $(this);
        evt.preventDefault();
        make_parent_active(self);

        iface._dataset = self.data('dataset');
        $('.page-header').text(iface._dataset.name);
        iface._render();
      });

      sample_container.find('a').each(function() {
        if(prev_dataset !== null && $(this).data('dataset').name === prev_dataset.name) {
          $(this).click();
        }
      });
    });
  });
}

google.load('visualization', '1.1', {packages: ['corechart', 'bar']});
google.setOnLoadCallback(main);

function main() {
  new Interface();
}
