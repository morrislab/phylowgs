function TreeUtil() {
}

TreeUtil.calc_radius = function(scale) {
  var min_area = 700, max_area = 8000;
  var area = Util.calc_in_range(min_area, max_area, scale);
  return Math.sqrt(area / Math.PI);
}
