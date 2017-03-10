function BestTreeFinder(tree_summ) {
  this._extract_indices(tree_summ);
}

BestTreeFinder.prototype._extract_indices = function(tree_summ) {
  this.linearity_indices = {};
  this.branching_indices = {};
  this.coclustering_indices = {};

  var self = this;
  Object.keys(tree_summ).forEach(function(tidx) {
    self.linearity_indices[tidx] = tree_summ[tidx].linearity_index;
    self.branching_indices[tidx] = tree_summ[tidx].branching_index;
    self.coclustering_indices[tidx] = tree_summ[tidx].coclustering_index;
  });
  this.num_trees = Object.keys(this.linearity_indices).length;
}

BestTreeFinder.prototype.calc_mean_index = function() {
  var mean = {};
  var self = this;
  ['linearity', 'branching', 'coclustering'].forEach(function(T) {
    var source = self[T + '_indices'];
    var sum = Object.keys(source).reduce(function(acc, tidx) {
      return acc + source[tidx];
    }, 0);
    mean[T] = sum / self.num_trees;
  });
  return mean;
}

BestTreeFinder.prototype._calc_euclid_dist = function(A, B) {
  var N = A.length;
  var dist = 0;
  for(var i = 0; i < N; i++) {
    dist += Math.pow(A[i] - B[i], 2);
  }
  return Math.sqrt(dist);
}

BestTreeFinder.prototype._calc_dists_from_median = function() {
  var N = this.num_trees;
  var dists = {};

  var self = this;
  Object.keys(this.linearity_indices).forEach(function(i) {
    var D = [];
    var vec_i = [self.linearity_indices[i], self.branching_indices[i], self.coclustering_indices[i]];

    Object.keys(self.linearity_indices).forEach(function(j) {
      var vec_j = [self.linearity_indices[j], self.branching_indices[j], self.coclustering_indices[j]];
      D.push(self._calc_euclid_dist(vec_i, vec_j));
    });

    var median = D.sort()[Math.floor(N/2)];
    dists[i] = median;
  });

  return dists;
}

BestTreeFinder.prototype._calc_dists_from_mean = function() {
  var dists = {};
  var mean_index = this.calc_mean_index();
  var meanvec = [mean_index.linearity, mean_index.branching, mean_index.coclustering];

  var self = this;
  Object.keys(this.linearity_indices).forEach(function(tidx) {
    var idxvec = [self.linearity_indices[tidx], self.branching_indices[tidx], self.coclustering_indices[tidx]];
    dists[tidx] = self._calc_euclid_dist(meanvec, idxvec);
  });
  return dists;
}

BestTreeFinder.prototype.calc_dists = function(normalize) {
  var dists = this._calc_dists_from_median();

  if(normalize) {
    // Since we use Euclidean distance and three-component vectors with elements // in [0, 1], the max distance is sqrt(3).
    var max_dist = Math.sqrt(3);
    var normed = [];
    Object.keys(dists).forEach(function(tidx) {
      normed[tidx] = dists[tidx] / max_dist;
    });
    dists = normed;
  }

  return dists;
}

BestTreeFinder.prototype.find_best_tree = function() {
  var min_dist = Number.POSITIVE_INFINITY;
  var best_tree_idx = null;
  var dists = this.calc_dists();

  Object.keys(dists).forEach(function(tidx) {
    var dist = dists[tidx];
    if(dist < min_dist) {
      min_dist = dist;
      best_tree_idx = parseInt(tidx, 10);
    }
  });

  return best_tree_idx;
}

