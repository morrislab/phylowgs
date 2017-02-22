Config = {
  font_size: 15
};

var MARKER = '<defs><marker id="head" orient="auto" markerWidth="2" markerHeight="4" refX="0" refY="2"><path d="M0,0 V4 L2,2 Z" fill="#000" /></marker></defs>';

function mean(list) {
  return list.reduce(function(a, b) { return a + b; }) / list.length;
}

function get_url_param(name) {
  name = name.replace(/[\[]/, "\\[").replace(/[\]]/, "\\]");
  var regex = new RegExp("[\\?&]" + name + "=([^&#]*)"),
      results = regex.exec(location.search);
  return results === null ? "" : decodeURIComponent(results[1].replace(/\+/g, " "));
}

function main() {
  new Interface();
  setTimeout(function() {
    StateManager.restore();
  }, 500);

  if(window.location.href.indexOf('debug=1') !== -1) {
    setTimeout(function() {
      $('#runs a:contains("cnvint")').click();
      $('#sample-list a:contains("0e7ac212")').click();
      $('.nav-tree-viewer').click();
    }, 1000);
  }
}

google.charts.load('45', {packages: ['corechart', 'bar', 'line']});
google.charts.setOnLoadCallback(main);
