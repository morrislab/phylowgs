google.load('visualization', '1.1', {packages: ['corechart', 'bar']});
google.setOnLoadCallback(main);

function draw_tree(root) {
  var m = [10, 20, 10, 20],
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

function display_tree(summary) {
  var adjlist = summary.structure;
  var pops = summary.populations;

  var max_area = 8000;
  var min_area = 700;
  var max_ssms = find_max_ssms(pops);

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
}

function make_parent_active(elem) {
  var parent = elem.parent();
  parent.siblings('li').removeClass('active');
  parent.addClass('active');
}

function sort_numeric(arr) {
  return arr.sort(function(a, b) {
    return parseInt(a, 10) - parseInt(b, 10);
  });
}

function mean(arr) {
  var sum = 0;
  for(var i = 0; i < arr.length; i++)
    sum += arr[i];
  return sum / arr.length;
}

function plot_vafs(muts) {
  var vafs = [];
  for(var ssm_id in muts.ssms) {
    var ssm = muts.ssms[ssm_id];
    var ssm_vafs = [];
    for(var i = 0; i < ssm.ref_reads.length; i++) {
      var a = ssm.ref_reads[i];
      var d = ssm.total_reads[i];
      ssm_vafs.push((d - a)/d);
    }
    vafs.push([mean(ssm_vafs)]);
  }

  var data = new google.visualization.DataTable();
  data.addColumn('number', 'VAF');
  data.addRows(vafs);

  var x_min = 0;
  var x_max = Math.max(1.0, array_max(vafs));
  var options = {
    title: 'VAFs (' + vafs.length + ' variants)',
    histogram: { bucketSize: 0.03 },
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

  var container = $('<div/>').appendTo('#container').get(0);
  var chart = new google.visualization.Histogram(container);
  chart.draw(data, options);
}

function array_max(arr) {
  return Math.max.apply(null, arr);
}

function render_tree(dataset) {
  var tree_container = $('#trees tbody');

  d3.json(dataset.summary_path, function(summary) {
    var tree_indices = sort_numeric(Object.keys(summary.trees));
    tree_container.empty();
    tree_indices.forEach(function(tidx) {
      var row = '<td class="tree-index">' + tidx + '</td><td class="tree-llh">' + summary.trees[tidx].llh.toFixed(1) + '</td>';
      $('<tr/>').html(row).appendTo(tree_container);
    });

    console.log($('#trees'));
    $('#trees').stupidtable();


    tree_container.find('tr').click(function(evt) {
      evt.preventDefault();
      var self = $(this);
      self.siblings().removeClass('active');
      self.addClass('active');

      var tidx = self.find('.tree-index').text();
      display_tree(summary.trees[tidx]);

      var summary_table = $('#tree-summary').show().find('tbody').empty();
      var pop_ids = sort_numeric(Object.keys(summary.trees[tidx].populations));
      pop_ids.forEach(function(pop_id) {
        var pop = summary.trees[tidx].populations[pop_id];
        var entries = [pop_id, pop.phi, pop.num_ssms, pop.num_cnvs].map(function(entry) {
          return '<td>' + entry + '</td>';
        });
        $('<tr/>').html(entries.join('')).appendTo(summary_table);
      });
    }).first().click();
    $('#tree-list').scrollTop(0);
  });
}

function extract_pops_with_top_phis(populations, desired_pops) {
  var pops = [];
  for(var pop_idx in populations) {
    pops.push(populations[pop_idx]);
  }
  pops.sort(function(a, b) {
    return parseInt(b.phi, 10) - parseInt(a.phi, 10);
  });

  // Exclude clonal population.
  var sliced = pops.slice(1, desired_pops + 1);
  while(sliced.length < desired_pops)
    sliced.push(null);
  return sliced;
}

function _render_phis(phis) {
  for(var i = 0; i < phis.length; i++) {
    var data = new google.visualization.DataTable();
    data.addColumn('number', 'Phi');
    data.addRows(phis[i]);

    var x_min = 0;
    var x_max = 1.0;
    var options = {
      title: 'Phi distribution (' + (i + 1) + ') (' + phis[i].length + ' values)',
      hAxis: {
        title: 'Phi',
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

function _render_ssm_counts(ssm_counts) {
  for(var i = 0; i < ssm_counts.length; i++) {
    var data = new google.visualization.DataTable();
    data.addColumn('number', 'SSMs');
    data.addRows(ssm_counts[i]);

    var options = {
      title: 'Number of SSMs (' + (i + 1) + ') (' + ssm_counts[i].length + ' values)',
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

function _render_vafs(dataset) {
  var muts_path = dataset.muts_path;
  d3.json(muts_path, function(muts) {
    plot_vafs(muts);
  });
}

function _render_pop_counts(pop_counts, min_ssms) {
  var histogram = {};
  var min_count = pop_counts.length, max_count = 0;
  console.log(pop_counts);
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
    hAxis: {
      title: 'Number of cancerous populations with at least one CNV or ' + min_ssms + ' SSMs',
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

function render_summary(dataset) {
  _render_vafs(dataset);

  var pops_to_examine = 3;
  var min_ssms = 3;

  var pop_counts = [];
  var phis = new Array(pops_to_examine);
  var ssm_counts = new Array(pops_to_examine);
  for(var i = 0; i < pops_to_examine; i++) {
    phis[i] = [];
    ssm_counts[i] = [];
  }

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

      var pops = extract_pops_with_top_phis(populations, pops_to_examine);
      for(var i = 0; i < pops.length; i++) {
        if(pops[i] !== null) {
          phis[i].push([pops[i].phi]);
          ssm_counts[i].push([pops[i].num_ssms]);
        }
      }
    }

    _render_phis(phis);
    _render_ssm_counts(ssm_counts);
    _render_pop_counts(pop_counts, min_ssms);
  });
}

function main() {
  var run_container = $('#runs');
  var sample_container = $('#samples');

  var renderer = render_tree;
  var dataset = null;

  $('.filter').keyup(function(evt) {
    var self = $(this);
    var filter_text = self.val();
    var elems = self.parents('.sidebar').find('.nav-sidebar').find('li');
    elems.hide();
    elems.filter(function() {
      return !($(this).text().indexOf(filter_text) === -1);
    }).show();
  });

  $('.navbar-nav a').click(function(evt) {
    evt.preventDefault();
    var self = $(this);
    make_parent_active(self);

    var mapping = {
      'nav-tree-summaries': render_summary,
      'nav-tree-viewer': render_tree
    };
    renderer = null;
    for(var nav_class in mapping) {
      if(self.hasClass(nav_class)) {
        renderer = mapping[nav_class];
        break;
      }
    }

    if(renderer !== null && dataset !== null) {
      $('#container').empty();
      $('#tree-summary').hide();
      renderer(dataset);
    }
  });

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

      dataset = null;
      sample_container.empty();

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


        dataset = self.data('dataset');
        if(renderer !== null && dataset !== null) {
          $('#container').empty();
          $('#tree-summary').hide();
          renderer(dataset);
        }
      });
    });
  });
}
