function TreePlotter() {
}

TreePlotter.prototype._sort_numeric = function(arr) {
  return arr.sort(function(a, b) {
    return parseInt(a, 10) - parseInt(b, 10);
  });
}

TreePlotter.prototype._calc_ccf = function(tree, pop_id) {
  var cp = tree.populations[pop_id].cellular_prevalence;
  if(parseInt(pop_id, 10) === 0)
    return cp.map(function(E) { return 0; });

  var ccf = [];
  // This only works for monoclonal trees, but it's the desired behaviour according to Quaid.
  var cellularity = tree.populations[1].cellular_prevalence;
  for(var i = 0; i < cellularity.length; i++) {
    ccf.push(cp[i] / cellularity[i]);
  }

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
      var pop_ids = tplotter._sort_numeric(Object.keys(summary.trees[tidx].populations));
      var root = tplotter._generate_tree_struct(summary.trees[tidx].structure, summary.trees[tidx].populations);
      tplotter._draw_tree(root);

      var summary_table = $('#snippets .tree-summary').clone().appendTo('#container');
      var num_pops = summary.trees[tidx].populations[pop_ids[0]].cellular_prevalence.length;
      summary_table.find('.cellprev').attr('colspan', num_pops);
      summary_table.find('.ccf').attr('colspan', num_pops);
      summary_table = summary_table.find('tbody');

      pop_ids.forEach(function(pop_id) {
        var pop = summary.trees[tidx].populations[pop_id];
        var cp = pop.cellular_prevalence.map(function(E) { return E.toFixed(3); });
        var ccf = tplotter._calc_ccf(summary.trees[tidx], pop_id);
        ccf = ccf.map(function(E) { return E.toFixed(3); });

        var entries = [pop_id].concat(cp).concat(ccf).concat([pop.num_ssms, pop.num_cnvs]).map(function(entry) {
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
