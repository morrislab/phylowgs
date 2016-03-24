function ClusterPlotter() {
  // horiz_padding should be set to the maximum radius of a node, so a node
  // drawn on a boundry won't go over the canvas edge. Since max_area = 8000,
  // we have horiz_padding = sqrt(8000 / pi) =~ 51.
  var horiz_padding = 51;
  this._M = [10, horiz_padding, 10, horiz_padding],
      this._width = 800 - this._M[1] - this._M[3],
      this._height = 600 - this._M[0] - this._M[2];

  $('#cluster-mode-chooser button').click(function() {
      var self = $(this);
      self.siblings('button').removeClass('active');
      self.addClass('active');
      // Trigger redraw.
      $('#clusters tr.active').click();
  });
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
    $('#cluster-size').stupidsort('desc');

    cluster_table.find('tr').click(function(evt) {
      evt.preventDefault();
      var self = $(this);
      self.siblings().removeClass('active');
      self.addClass('active');

      var cidx = self.find('.cluster-index').text();

      var active_cluster_mode = $('#cluster-mode-chooser .active');
      if(active_cluster_mode.hasClass('centroid-chooser')) {
        cplotter._show_centroid(formatted_clusters[cidx]);
      } else if(active_cluster_mode.hasClass('average-chooser')) {
        cplotter._show_average(formatted_clusters[cidx]);
      } else {
        throw 'No cluster mode chosen';
      }
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

ClusterPlotter.prototype._calc_max_spanning_tree = function(cluster) {
  // Use Kruskal's algorithm to get max spanning tree.
  var edges = [];
  var vertex_sets = [];
  var mst = {};

  Object.keys(cluster.populations).forEach(function(popidx) {
    var popidx = parseInt(popidx, 10);
    var pop = cluster.populations[popidx];
    vertex_sets.push([popidx]);

    Object.keys(pop.children).forEach(function(childidx) {
      var childidx = parseInt(childidx, 10);
      var weight = pop.children[childidx];
      edges.push({weight: weight, from: popidx, to: childidx});
    });
  });

  edges.sort(function(a, b) {
    // Sort in descending order.
    return b.weight - a.weight;
  });

  var _get_set_idx = function(vert) {
    for(var i = 0; i < vertex_sets.length; i++) {
      if(typeof vertex_sets[i] !== 'undefined' && vertex_sets[i].indexOf(vert) > -1) {
        return i;
      }
    }
    return -1;
  };

  edges.forEach(function(edge) {
    var from_set_idx = _get_set_idx(edge.from);
    var to_set_idx = _get_set_idx(edge.to);
    if(from_set_idx === to_set_idx) {
      return;
    }

    vertex_sets[from_set_idx] = vertex_sets[from_set_idx].concat(vertex_sets[to_set_idx]);
    delete vertex_sets[to_set_idx];

    if(!mst.hasOwnProperty(edge.from)) {
      mst[edge.from] = [];
    }
    mst[edge.from].push(edge.to);
  });

  return mst;
}

ClusterPlotter.prototype._create_pops = function(cluster) {
  var pops = {};

  Object.keys(cluster.populations).forEach(function(popidx) {
    var pop = cluster.populations[popidx];
    var cell_prevs = Util.transpose(pop.cellular_prevalences).map(function(cp_estimates) {
      return Util.mean(cp_estimates);
    });
    pops[popidx] = {
      num_ssms: Math.round(Util.mean(pop.num_ssms)),
      // TODO: actually provide CNV counts rather than dummy 0 values.
      num_cnvs: 0,
      cellular_prevalence: cell_prevs,
    };
  });
  return pops;
}

ClusterPlotter.prototype._show_centroid = function(cluster) {
  var exemplar_structure = this._calc_max_spanning_tree(cluster);
  var exemplar_pops = this._create_pops(cluster);

  var tplotter = new TreePlotter();
  tplotter.draw(exemplar_pops, exemplar_structure);
}

ClusterPlotter.prototype._show_average = function(cluster) {
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
