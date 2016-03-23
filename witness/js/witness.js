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

google.load('visualization', '1.1', {packages: ['corechart', 'bar']});
google.setOnLoadCallback(main);

function main() {
  new Interface();
  /*setTimeout(function() {
    $('#runs a:contains("steph")').click();
    $('#sample-list a:contains("SJBALL022609")').click();
    $('.nav-clustered-trees').click();
  }, 1000);*/
}
