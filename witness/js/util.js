function Util() {
}

Util.mean = function(list) {
  return list.reduce(function(a, b) { return a + b; }) / list.length;
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
