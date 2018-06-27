function SSM_Viewer(dataset){
  this.muts = null;
  this.ssm_ids = null;
  this.ssm_names = null;
  this.num_filled_rows = 0;
  this.dataset = dataset;

  var self = this;
  d3.json(dataset.muts_path, function(muts){
    self.muts = muts;
  })
}

SSM_Viewer.prototype.render = function(container, tidx, num_samples, sample_names){
  var ssm_viewer_table = $('#ssm-viewer-table').clone().appendTo(container);
  var ssm_info_row = ssm_viewer_table.find('.ssm-info-row');

  this._setup_table(ssm_viewer_table, ssm_info_row, num_samples, sample_names);
  this._setup_ssm_entry_event(ssm_info_row, tidx, num_samples);

  if(this.ssm_ids !== null){
    this._update_row(ssm_info_row, tidx, num_samples);
  }
}

SSM_Viewer.prototype._setup_table = function(table, ssm_info_row, num_samples, sample_names){
  // SET UP THE HEADER
  var header = table.find('thead');
  // Empty cells for SSM, ID, Name and Node columns
  var samps_header = ['&mdash;', '&mdash;', '&mdash;', '&mdash;'];
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
    $('<th class=vaf-samp' + i + ' title=""> </th>').html('').appendTo(ssm_info_row);
}

SSM_Viewer.prototype._setup_ssm_entry_event = function(ssms_info_row, tidx, num_samples){
  //What to do once enter is hit after filling in the ssm textbox
  var sViewer = this;
  var input = ssms_info_row.find('.ssm-input');
  input.keypress(function(evt) {
    var key = evt.which;
    if(key==13){
      var self = $(this);
      var ssm_text = self.val();
      //Determine the ID and the name of the ssm from the input. Could have input either the ID itself or the name of the ssm
      if(ssm_text in sViewer.muts.ssms){
        sViewer.ssm_ids = ssm_text;
        if(sViewer.muts.ssms[ssm_text].hasOwnProperty('name'))
          sViewer.ssm_names = sViewer.muts.ssms[ssm_text].name;
        else
          //insert a dash if the name isn't available
          sViewer.ssm_names = '\u2014'; 
      }else{
        Object.keys(sViewer.muts.ssms).forEach(function(sidx){
          if(sViewer.muts.ssms[sidx].hasOwnProperty("name") & (sViewer.muts.ssms[sidx].name == ssm_text)){
            sViewer.ssm_ids = sidx;
            sViewer.ssm_names = ssm_text;
            return;
          }
        })
      }
      sViewer._update_row(ssms_info_row, tidx, num_samples);
    }
  });
}

SSM_Viewer.prototype._update_row = function(ssms_info_row, tidx, num_samples){
  var sViewer = this;
  //Clear the input box.
  ssms_info_row.find('.ssm-input').val('');
  //Set the ID and Name columns
  ssms_info_row.find('.ssms-id').text(this.ssm_ids);
  ssms_info_row.find('.ssms-name').text(this.ssm_names);
  //Load in the mutass json, find which subpopulation the ssm belongs to, set the table text
  var mutass_treepath = sViewer.dataset.mutass_path + "/" + tidx  + '.json';
  d3.json(mutass_treepath, function(mutass){
    var pop = sViewer._find_pop_containing_ssm(mutass, sViewer.ssm_ids);
    ssms_info_row.find('.ssms-node').text(pop);
  })
  //Get the rest of the info from the muts file and fill in the row
  if ((typeof this.muts.ssms[this.ssm_ids]) !== 'undefined'){
    var tot_reads = this.muts.ssms[this.ssm_ids].total_reads;
    var ref_reads = this.muts.ssms[this.ssm_ids].ref_reads;
    for(i=0;i<num_samples;i++){
      var vaf = (tot_reads[i] - ref_reads[i]) / tot_reads[i];
      var var_reads = tot_reads[i] - ref_reads[i];
      //table entry text
      ssms_info_row.find('.vaf-samp' + i).text(vaf.toFixed(3));
      //tooltip text
      ssms_info_row.find('.vaf-samp' + i).attr("title","Variant Reads=" + var_reads + "; Total Reads=" + tot_reads[i]);
    }
  }else{
    for(i=0;i<num_samples;i++)
      ssms_info_row.find('.vaf-samp' + i).text('');
      ssms_info_row.find('.vaf-samp' + i).attr("title",'');
  }
}

SSM_Viewer.prototype._find_pop_containing_ssm = function(mutass,ssm_ids){
  var ssms_pop = -1;
  //Search through the json for the ssm, report which pop it's in. If none, do nothing.
  Object.keys(mutass.mut_assignments).forEach(function(pop){
    if(ssms_pop !== -1)
      return;
    var this_pop = mutass.mut_assignments[pop];
    Object.keys(this_pop.ssms).forEach(function(ssm){
      //Check the key and, if it exists, the name of the ssm.
      if(this_pop.ssms[ssm] == ssm_ids){
        ssms_pop = pop;
        return;
      }
    })
  })
  return ssms_pop;
}