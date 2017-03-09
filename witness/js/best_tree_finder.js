function BestTreeFinder(lin_idx, branch_idx, tree_idx) {
  this._mean_index = null;
  this._lin_idx = lin_idx;
  this._branch_idx = branch_idx;
  this._tree_idx = tree_idx;
}

BestTreeFinder.prototype.calc_mean_index = function() {
  if(this._mean_index === null)
    this.make_indices();
  return this._mean_index;
}

BestTreeFinder.prototype.make_indices = function() {
  var N = this._lin_idx.length;
  var indices = [];
  var mean_indices = [0, 0, 0];

  for(var i = 0; i < N; i++) {
    var cocluster_idx = 1 - (this._lin_idx[i] + this._branch_idx[i]);
    indices.push([this._lin_idx[i], this._branch_idx[i], cocluster_idx]);
    mean_indices[0] += this._lin_idx[i];
    mean_indices[1] += this._branch_idx[i];
    mean_indices[2] += cocluster_idx;
  }
  mean_indices[0] /= N;
  mean_indices[1] /= N;
  mean_indices[2] /= N;
  this._mean_index = mean_indices;

  return indices;
}

BestTreeFinder.prototype._calc_euclid_dist = function(A, B) {
  var N = A.length;
  var dist = 0;
  for(var i = 0; i < N; i++) {
    dist += Math.pow(A[i] - B[i], 2);
  }
  return Math.sqrt(dist);
}

BestTreeFinder.prototype.calc_dists_from_mean = function(normalize) {
  var dists = {};
  var indices = this.make_indices();

  var K = this._mean_index.length;
  var max_len = Math.sqrt(K);

  var self = this;
  indices.forEach(function(idxs, i) {
    dists[i] = self._calc_euclid_dist(idxs, self._mean_index);
    // Don't normalize by default -- caller must explicitly request this.
    if(normalize) {
      dists[i] /= max_len;
    }
    if(parseInt(self._tree_idx[i], 10) !== i) {
      throw "tree_idx doesn't match expected index " + parseInt(self._tree_idx[i], 10) + ", " + i;
    }
  });
  return dists;
}

BestTreeFinder.prototype.find_best_tree = function() {
  var self = this;
  var min_dist = Number.POSITIVE_INFINITY;
  var best_tree_idx = null;
  var indices = this.make_indices();
  var dists = this.calc_dists_from_mean();

  Object.keys(dists).forEach(function(tidx) {
    var dist = dists[tidx];
    if(dist < min_dist) {
      min_dist = dist;
      best_tree_idx = parseInt(tidx, 10);
    }
  });
  return best_tree_idx;
}

