function TreeViewer() {
}

TreeViewer.prototype._plot_pop_vafs = function(dataset, tidx) {
  if(!(dataset.hasOwnProperty('muts_path') && dataset.hasOwnProperty('mutass_path')))
    return;

  var pop_vaf_plotter = new PopVafPlotter();
  d3.json(dataset.muts_path, function(muts) {
    d3.json(dataset.mutass_path + '/' + tidx + '.json', function(mutass) {
      pop_vaf_plotter.plot(muts, mutass);
    });
  });
}

TreeViewer.prototype.render = function(dataset) {
  $('#tree-list').show();
  var tree_container = $('#trees tbody');

  var tplotter = this;
  d3.json(dataset.summary_path, function(summary) {
    var tree_indices = Util.sort_ints(Object.keys(summary.trees));
    tree_container.empty();

    var first_tree_idx = tree_indices[0];
    var first_pop_idx = Object.keys(summary.trees[first_tree_idx].populations)[0];
    var num_samples = summary.trees[first_tree_idx].populations[first_pop_idx].cellular_prevalence.length;

    tree_indices.forEach(function(tidx) {
      var total_ssms = 0;
      Object.keys(summary.trees[tidx].populations).forEach(function(pidx) {
        total_ssms += summary.trees[tidx].populations[pidx].num_ssms;
      });

      var normllh_nats = -summary.trees[tidx].llh / total_ssms;
      normllh_nats /= num_samples;
      var normllh_bits = normllh_nats / Math.log(2);

      var row = '<td class="tree-index">' + tidx + '</td>'
        + '<td class="tree-llh">' + normllh_bits.toFixed(1) + '</td>'
        + '<td class="tree-nodes">' + Object.keys(summary.trees[tidx].populations).length + '</td>';
      ['linearity_index', 'branching_index', 'clustering_index'].forEach(function(idxname) {
        var val = summary.trees[tidx].hasOwnProperty(idxname) ? summary.trees[tidx][idxname].toFixed(2) : '&mdash;';
        row += '<td>' + val + '</td>';
      });
      $('<tr/>').html(row).appendTo(tree_container);
    });

    $('#trees').stupidtable();

    var already_autosorted = false;
    $('#trees').bind('aftertablesort', function() {
      if(already_autosorted)
        return;
      // If any restore events are in progress (e.g., waiting on a timer),
      // don't automatically click on the first element.
      if(StateManager.restoring > 0)
        return;
      tree_container.find('tr:first').click();
      already_autosorted = true;
    });

    // If direction not specified, this can end up being ascending or
    // descending sort, depending on prior sort state of table.
    $('#tree-llh').stupidsort('asc');

    tree_container.find('tr').click(function(evt) {
      evt.preventDefault();
      var self = $(this);
      self.siblings().removeClass('active');
      self.addClass('active');

      var tidx = self.find('.tree-index').text();
      StateManager.update('tidx', tidx);
      var tree_plotter = new TreePlotter();

      tree_plotter.draw(summary.trees[tidx].populations, summary.trees[tidx].structure, summary.trees[tidx].root, summary.params);
      tplotter._plot_pop_vafs(dataset, tidx);
    });
    $('#tree-list').scrollTop(0);
  });
}
