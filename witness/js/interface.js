function make_parent_active(elem) {
  var parent = elem.parent();
  parent.siblings('li').removeClass('active');
  parent.addClass('active');
}

function Interface() {
  this._activate_filters();
  this._activate_navbar();
  this._load_samples();

  this._available_renderers = {
    'summarizer': new TreeSummarizer(),
    'tree_viewer': new TreeViewer(),
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
      'nav-tree-viewer': 'tree_viewer',
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
