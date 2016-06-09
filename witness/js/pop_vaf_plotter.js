function PopVafPlotter() {
}

PopVafPlotter.prototype._compile_mut_vafs = function(muts, mutass) {
  var mut_vafs = {};

  Object.keys(mutass).forEach(function(popidx) {
    var ssms = mutass[popidx]['ssms'];
    if(ssms.length === 0)
      return;
    mut_vafs[popidx] = ssms.map(function(ssm_id) {
      var ref_reads = muts.ssms[ssm_id].ref_reads;
      var total_reads = muts.ssms[ssm_id].total_reads;
      var vaf = [];
      for(var i = 0; i < ref_reads.length; i++) {
        vaf.push((total_reads[i] - ref_reads[i]) / total_reads[i]);
      }
      return vaf;
    });
  });

  return mut_vafs;
}

PopVafPlotter.prototype._make_traces = function(mut_vafs, popidx) {
  var first_pop_idx = Object.keys(mut_vafs)[0]; // Unsorted, but doesn't matter
  var num_samples = mut_vafs[first_pop_idx][0].length;

  var traces = [];
  for(var traceidx = 0; traceidx < num_samples; traceidx++) {
    traces.push(mut_vafs[popidx].map(function(vafs_for_pop) {
      if(vafs_for_pop.length !== num_samples)
        throw 'Expected ' + num_samples + ' values, but got ' + vafs_for_pop.length;
      return vafs_for_pop[traceidx];
    }));
  }
  return traces;
}

PopVafPlotter.prototype._plot_histogram = function(popidx, trace_data) {
  var traces = trace_data.map(function(data, idx) {
    return {
      x: data,
      name: 'Sample ' + (idx + 1),
      type: 'histogram',
      opacity: 0.75
    };
  });

  var layout = {
    title: 'VAFs for population ' + popidx,
    barmode: 'overlay',
    xaxis: { title: 'VAFs'},
    yaxis: { title: 'Number of SSMs'}
  };

  var container = document.querySelector('#container');
  var plot_container = document.createElement('div');
  container.appendChild(plot_container);
  Plotly.newPlot(plot_container, traces, layout);
}

PopVafPlotter.prototype.plot = function(muts, mutass) {
  var mut_vafs = this._compile_mut_vafs(muts, mutass['mut_assignments']);

  var popidxs = Object.keys(mut_vafs).sort(function(a, b) {
    return parseInt(a, 10) - parseInt(b, 10);
  });
  var self = this;
  popidxs.forEach(function(popidx) {
    var trace_data = self._make_traces(mut_vafs, popidx);
    self._plot_histogram(popidx, trace_data);
  });
}

