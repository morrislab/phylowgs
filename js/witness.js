function draw_tree(root) {
  var m = [10, 20, 10, 20],
      w = 800 - m[1] - m[3],
      h = 400 - m[0] - m[2],
      i = 0;

  var tree = d3.layout.tree()
      .size([h, w])
      .sort(function(a, b) {
        return d3.ascending(parseInt(a.name, 10), parseInt(b.name, 10));
      });

  var diagonal = d3.svg.diagonal()
      .projection(function(d) { return [d.y, d.x]; });

  var vis = d3.select('#container').append('svg:svg')
      .attr('width', w + m[1] + m[3])
      .attr('height', h + m[0] + m[2])
      .append('svg:g')
      .attr('transform', 'translate(' + m[3] + ',' + m[0] + ')');

  // Compute the new tree layout.
  var nodes = tree.nodes(root);

  // Normalize for fixed-depth.
  nodes.forEach(function(d) { console.log([d.depth, d.y / d.depth]); d.y = d.depth * 180; });

  // Update the nodes…
  var node = vis.selectAll('g.node')
      .data(nodes, function(d) { return d.name; });

  // Enter any new nodes at the parent's previous position.
  var nodeEnter = node.enter().append('svg:g')
      .attr('class', 'node');

  nodeEnter.append('svg:circle')
      .attr('r', function(d) { return d.radius; })
      .style('fill', 'lightsteelblue');

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
      .data(tree.links(nodes), function(d) { return d.target.name; });

  // Enter any new links at the parent's previous position.
  link.enter().insert('svg:path', 'g')
      .attr('class', 'link');

  // Transition links to their new position.
  link.attr('d', diagonal);

  // Transition exiting nodes to the parent's new position.
  link.exit().remove();
}

function get_url_param(name) {
  name = name.replace(/[\[]/, "\\[").replace(/[\]]/, "\\]");
  var regex = new RegExp("[\\?&]" + name + "=([^&#]*)"),
      results = regex.exec(location.search);
  return results === null ? "" : decodeURIComponent(results[1].replace(/\+/g, " "));
}

function find_max_ssms(populations) {
  var max_ssms = 0;
  for(var pop_id in populations) {
    var pop = populations[pop_id];
    if(pop.num_ssms > max_ssms)
      max_ssms = pop.num_ssms;
  }
  return max_ssms;
}

function main() {
  var data_path = 'data/' + get_url_param('summary') + '.summ.json.gz';
  d3.json(data_path, function(summary) {
    var tree_id = 2001;
    var adjlist = summary.trees[tree_id].structure;
    var pops = summary.trees[tree_id].populations;

    var max_area = 8000;
    var min_area = 700;
    var max_ssms = find_max_ssms(summary.trees[tree_id].populations);

    var _add_node = function(node_id, struct) {
      struct.name = node_id;

      var num_ssms = pops[node_id]['num_ssms'];
      var area = min_area + (num_ssms / max_ssms)*(max_area - min_area);
      if(area < min_area) area = min_area;
      if(area > max_area) area = max_area;
      struct.radius = Math.sqrt(area / Math.PI);

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
    _add_node(0, root);
    draw_tree(root);
  });
}

main()
