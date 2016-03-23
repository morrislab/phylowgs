function TreeViewer() {
}


TreeViewer.prototype.render = function(dataset) {
  $('#tree-list').show();
  var tree_container = $('#trees tbody');

  var tplotter = this;
  d3.json(dataset.summary_path, function(summary) {
    var tree_indices = Util.sort_ints(Object.keys(summary.trees));
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
      var tree_plotter = new TreePlotter();
      tree_plotter.draw(summary.trees[tidx].populations, summary.trees[tidx].structure);
    });
    $('#tree-list').scrollTop(0);
  });
}
