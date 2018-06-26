function SSM_Viewer(dataset){
  this.muts = null;
  this.ssm_id = null;
  this.dataset = dataset;

  if(this.muts === null){
    self = this;
    d3.json(dataset.muts_path, function(muts){
      self.muts = muts;
    })
  }
}

SSM_Viewer.prototype._find_pop_with_ssm = function(mutass,ssm_id){
  var ssms_pop = -1;
  //Search through the json for the ssm, report which pop it's in. If none, do nothing.
  Object.keys(mutass.mut_assignments).forEach(function(pop){
    if(ssms_pop !== -1)
      return;
    var this_pop = mutass.mut_assignments[pop];
    Object.keys(this_pop.ssms).forEach(function(ssm){
      //Check the key and, if it exists, the name of the ssm.
      if(this_pop.ssms[ssm] == ssm_id){
        ssms_pop = pop;
        return;
      }
    })
  })
  return ssms_pop;
}

SSM_Viewer.prototype._update_table = function(ssms_info_row, tidx, num_samples){
  var ssms_node_txtbox = ssms_info_row.find('.ssms-node');
  var sViewer = this;
  //Make sure the ssm id is in the input box
  var ssm_input = ssms_info_row.find('.ssm-input');
  ssm_input.val(sViewer.ssm_id);
  //Load in the mutass json, find which subpopulation the ssm belongs to, set the table text
  var mutass_treepath = sViewer.dataset.mutass_path + "/" + tidx  + '.json';
  d3.json(mutass_treepath, function(mutass){
    var pop = sViewer._find_pop_with_ssm(mutass, sViewer.ssm_id);
    ssms_node_txtbox.text(pop);
  })
  //Get the rest of the info from the muts file and fill in the row
  if ((typeof this.muts.ssms[this.ssm_id]) !== 'undefined'){
    var tot_reads = this.muts.ssms[this.ssm_id].total_reads;
    var ref_reads = this.muts.ssms[this.ssm_id].ref_reads;
    for(i=0;i<num_samples;i++){
      var vaf = (tot_reads[i] - ref_reads[i]) / tot_reads[i];
      var var_reads = tot_reads[i] - ref_reads[i];
      //table entry text
      ssms_info_row.find('#vaf-samp' + i).text(vaf.toFixed(3));
      //Tooltip text
      ssms_info_row.find('#vaf-samp' + i).attr("title","Variant Reads=" + var_reads + "; Total Reads=" + tot_reads[i]);
    }
  }else{
    for(i=0;i<num_samples;i++)
      ssms_info_row.find('#vaf-samp' + i).text('');
      ssms_info_row.find('#vaf-samp' + i).attr("title",'');
  }
}

SSM_Viewer.prototype._setup_ssm_entry_event = function(input, tidx, ssms_info_row, num_samples){
  //What to do once enter is hit after filling in the ssm textbox
  sviewer = this;
  input.keypress(function(evt) {
    var key = evt.which;
    if(key==13){
      var self = $(this);
      var ssm_text = self.val();
      //Determine the ID of the ssm from the input. Could have input either the ID itself or the name of the ssm
      if(ssm_text in sviewer.muts.ssms){
        sviewer.ssm_id = ssm_text;
      }else{
        Object.keys(sviewer.muts.ssms).forEach(function(sidx){
          if(sviewer.muts.ssms[sidx].hasOwnProperty("name") & (sviewer.muts.ssms[sidx].name == ssm_text)){
            sviewer.ssm_id = sidx;
            return;
          }
        })
      }
      sviewer._update_table(ssms_info_row, tidx, num_samples);
    }
  });
}

SSM_Viewer.prototype._setup_table = function(table, num_samples, sample_names, ssm_info_row){

  // SET UP THE HEADER
  var header = table.find('thead');
  // Empty cells for SSM and node columns
  var samps_header = ['&mdash;', '&mdash;'];
  // Sample names for the rest of the columns
  samps_header = samps_header.concat(sample_names);
  var samps_header = samps_header.map(function(entry) {
    return '<th>' + entry + '</th>';
  });
  table.find('.vafs-header').attr('colspan', num_samples);
  $('<tr/>').html(samps_header.join('')).appendTo(header);

  // SET UP THE BODY
  // Add columns with necessary ids for each sample in Ref Reads and Tot Reads
  for(i=0;i<num_samples;i++)
    //note:the title attribute can be used as a tooltip
    $('<th id=vaf-samp' + i + ' title=""> </th>').html('').appendTo(ssm_info_row);
}

SSM_Viewer.prototype.render = function(container, num_samples, sample_names, tidx){
  var ssm_viewer_table = $('#ssm-viewer-table').clone().appendTo(container);
  var ssm_info_row = ssm_viewer_table.find('.ssm-info');
  var ssm_input = ssm_viewer_table.find('.ssm-input');

  this._setup_table(ssm_viewer_table, num_samples, sample_names, ssm_info_row);
  this._setup_ssm_entry_event(ssm_input, tidx, ssm_info_row, num_samples);

  if(this.ssm_id !==null){
    this._update_table(ssm_info_row, tidx, num_samples);
  }
}