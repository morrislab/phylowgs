function SSM_Viewer(){
}

SSM_Viewer.prototype._find_pop_with_ssm = function(mutass,ssm_id){
  ssms_pop = -1
  //Search through the json for the ssm, report which pop it's in. If none, do nothing.
  Object.keys(mutass.mut_assignments).forEach(function(pop){
    if(ssms_pop !== -1)
      return
    this_pop = mutass.mut_assignments[pop];
    Object.keys(this_pop.ssms).forEach(function(ssm){
      //console.log(this_pop.ssms[ssm])
      if(this_pop.ssms[ssm] == ssm_id){
        ssms_pop = pop;
        return;
      }
    })
  })
  return ssms_pop;
}

SSM_Viewer.prototype._setup_ssm_entry_event = function(input, tidx, dataset, ssms_info_row, num_samples){
  var ssms_node_txtbox = ssms_info_row.find('.ssms-node');
  var sViewer = this;
  //What to do once enter is hit after filling in the ssm textbox
  input.keypress(function(evt) {  
    var key = evt.which;
    if(key==13){
      var self = $(this);
      var ssm_text = self.val();
      //Load in the mutass json, find which subpop the ssm belongs to, set the table text
      var mutass_treepath = dataset.mutass_path + "/" + tidx  + '.json'
      d3.json(mutass_treepath, function(mutass){
        var pop = sViewer._find_pop_with_ssm(mutass, ssm_text);
        ssms_node_txtbox.text(pop);
      })
      //Get the rest of the info from the muts file and fill in the row
      d3.json(dataset.muts_path, function(muts){
        if ((typeof muts.ssms[ssm_text]) !== 'undefined'){
          var tot_reads = muts.ssms[ssm_text].total_reads;
          var ref_reads = muts.ssms[ssm_text].ref_reads;
          for(i=0;i<num_samples;i++){
            ssms_info_row.find('.ref-read' + i).text(ref_reads[i]);
            ssms_info_row.find('.tot-read' + i).text(tot_reads[i]);
          }
        }else{
          for(i=0;i<num_samples;i++){
            ssms_info_row.find('.ref-read' + i).text('');
            ssms_info_row.find('.tot-read' + i).text('');
          }
        }
        
      })
    }
  });
}

SSM_Viewer.prototype._setup_table = function(table, num_samples, sample_names, ssm_info_row){
  table.find('.refreads').attr('colspan', num_samples);
  table.find('.totreads').attr('colspan', num_samples);
  var header = table.find('thead');
  var body = table.find('tbody');

  // SETUP THE HEADER
  // Empty cells for SSM, and node columns.
  var samps_header = ['&mdash;', '&mdash;']
  // Sample names for the rest of the columns
  samps_header = samps_header.concat(sample_names).concat(sample_names);
  var samps_header = samps_header.map(function(entry) {
    return '<th>' + entry + '</th>';
  });
  $('<tr/>').html(samps_header.join('')).appendTo(header);

  // SETUP THE BODY
  // Add columns with necessary ids for each sample in Ref Reads and Tot Reads
  for(i=0;i<num_samples;i++)
    $('<th class=ref-read' + i + '></th>').html('').appendTo(ssm_info_row);
  for(i=0;i<num_samples;i++)
    $('<th class=tot-read' + i + '></th>').html('').appendTo(ssm_info_row);
}

SSM_Viewer.prototype.render = function(container, num_samples, sample_names, tidx, dataset){
  /*TODO:
  - Finish setting up the table.
    - (done) Automatically expand the table columns to account for the number of samples.
    - (done) Can type in ssm name.
    - Set up the enter event for the ssm-input box, have it fill in the columns of the table.
    - (Optional) Add more rows to the table? That way can view multiple ssms at the same time.
  - Load in the mutass file corresponding to the selected tree. New server handler should be able to do it,
    just have to send the correct signal to fetch the file, which im not sure how to do right now
  - Will also need access to the muts file for the reference reads and varient reads. Can probably
    just look at other parts of witness code to see where muts is used.
  - 
  */
  
  var ssm_viewer_table = $('#ssm-viewer-table').clone().appendTo(container);
  var ssm_info_row = ssm_viewer_table.find('.ssm-info')
  var ssm_input = ssm_viewer_table.find('.ssm-input');

  this._setup_table(ssm_viewer_table, num_samples, sample_names, ssm_info_row);
  this._setup_ssm_entry_event(ssm_input, tidx, dataset, ssm_info_row, num_samples);
  
}