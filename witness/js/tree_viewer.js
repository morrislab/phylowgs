function TreeViewer() {
}

TreeViewer.prototype._plot_pop_vafs = function(dataset, tidx) {
  if(!(dataset.hasOwnProperty('muts_path') && dataset.hasOwnProperty('mutass_path')))
    return;

  var pop_vaf_plotter = new PopVafPlotter();
  d3.json(dataset.muts_path, function(muts) {
    d3.json(dataset.mutass_path + '/' + tidx + '.json', function(mutass) {
      pop_vaf_plotter.plot(muts, mutass);
    })
  })
}

TreeViewer.prototype.render = function(dataset) {
  var tviewer = this;
  d3.json(dataset.summary_path, function(summary) {
    //Fill in tree table.
    $('#tree-list').show();
    var tree_container = $('#trees tbody');
    tree_container.empty();
    var tree_table = new TreeTable();
    tree_table.create(summary, tree_container);
    $('#trees').stupidtable();
    
    //Display the info for the first tree
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
    })

    // If direction not specified, this can end up being ascending or
    // descending sort, depending on prior sort state of table.
    $('#tree-llh').stupidsort('asc');

    //Define click event for rows in tree table.
    tree_container.find('tr').click(function(evt) {
      evt.preventDefault();
      var self = $(this);
      self.siblings().removeClass('active');
      self.addClass('active');
      
      var tidx = self.find('.tree-index').text();
      StateManager.update('tidx', tidx);
      var tree_plotter = new TreePlotter();
      var container = document.querySelector('#container');

      tree_plotter.draw(tidx, summary.trees[tidx].populations, summary.trees[tidx].structure, summary.trees[tidx].root, summary.params, dataset);
      tviewer._plot_pop_vafs(dataset, tidx);
    })
    $('#tree-list').scrollTop(0);
  })
}
