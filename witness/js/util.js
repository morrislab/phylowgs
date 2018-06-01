function Util() {
}

Util.mean = function(list) {
  return list.reduce(function(a, b) { return a + b; }) / list.length;
}

Util.stdev = function(list) {
  var mean = Util.mean(list);
  var diffs = list.map(function(E) { return Math.pow(E - mean, 2); });
  return Math.sqrt(Util.mean(diffs));
}

Util.array_max = function(arr) {
  return Math.max.apply(null, arr);
}

// scale must be in [0, 1].
Util.calc_in_range = function(min, max, scale) {
  return min + scale*(max - min);
}

Util.sort_ints = function(arr) {
  return arr.sort(function(a, b) {
    return parseInt(a, 10) - parseInt(b, 10);
  });
}

Util.transpose = function(arr) {
  // Taken from http://stackoverflow.com/a/17428705
  return arr[0].map(function(_, i) {
    return arr.map(function(row) {
      return row[i];
    });
  });
}

Util.have_cluster_data = function(summary){
  // Check to see if using older data where no clustering was done. If so, may not be able to plot some graphs.
  var have_clust_data = false;
  console.log(summary)
  if(summary.hasOwnProperty('clusters')){
    if(summary.clusters[Object.keys(summary.clusters)[0]].hasOwnProperty('members')){
      var have_clust_data = true;
    }
  };
  //Object.keys(summary.clusters).forEach(function(cKey){
  //  if (!summary.clusters[cKey].hasOwnProperty('members')){
  //    using_old_data = true;
  //    return
  //  }
  //})
  return have_clust_data;
}
