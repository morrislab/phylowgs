function SSM_Viewer(dataset){
  //TODO:
  // - just calc num_samples from sample_names
  // - Leave original input in place, hide any others from future rows.
  // - When new information is input, just add a new row (unless ssm already exists?!?)
  // - 
  this._muts = null;
  this._ssm_ids = [];
  this._ssm_names = [];
  this._dataset = dataset;
  this._table = null;
  this._input = null;
  this._first_row = null;

  var self = this;
  d3.json(dataset.muts_path, function(muts){
    self._muts = muts;
  })

  this.render = function(container, tidx, num_samples, sample_names){
    this._table = $('#ssm-viewer-table').clone().appendTo(container);
    this._first_row = this._table.find('.ssm-info-row');
    this._input = this._first_row.find('.ssm-input');

    this._setup_table(num_samples, sample_names);
    this._setup_ssm_entry_event(tidx, num_samples);

    if(this._ssm_ids.length > 0){
      this._update_table(tidx, num_samples);
    }
  }

  this._setup_table = function(num_samples, sample_names){
    //Make use of the table as defined in index.html, add additional columns 
    //according to the number of samples in this dataset. Add ids as necessary
    //to allow straight forward insertion of ssm data.

    // SET UP THE HEADER
    var header = this._table.find('thead');
    // Empty cells for SSM, ID, Name and Node columns
    var samps_header = ['&mdash;', '&mdash;', '&mdash;', '&mdash;'];
    // Sample names for the rest of the columns
    samps_header = samps_header.concat(sample_names);
    var samps_header = samps_header.map(function(entry) {
      return '<th>' + entry + '</th>';
    });
    this._table.find('.vafs-header').attr('colspan', num_samples);
    $('<tr/>').html(samps_header.join('')).appendTo(header);

    // SET UP THE FIRST ROW OF THE BODY
    // Add columns with necessary ids for each sample in Ref Reads and Tot Reads
    for(i=0;i<num_samples;i++)
      //note:the title attribute can be used as a tooltip
      $('<th id=vaf-samp' + i + '-0 title=""> </th>').html('').appendTo(this._first_row);
  }

  this._setup_ssm_entry_event = function(tidx, num_samples){
    //What to do once enter is hit after filling in the ssm input box
    var sViewer = this;
    this._input.keypress(function(evt) {
      var key = evt.which;
      if(key==13){
        var self = $(this);
        var ssm_text = self.val();
        //Determine the ID and the name of the ssm from the input. Could have input either the ID itself or the name of the ssm
        var correct_input = sViewer._set_ID_and_name_from_input_text(ssm_text);
        if(correct_input){
          sViewer._add_row(num_samples, sViewer._ssm_ids.length-1);
          sViewer._fill_row(tidx, num_samples, sViewer._ssm_ids.length-1);
        }
      }
    });
  }

  this._update_table = function(tidx,num_samples){
    //If there has already been ssm IDs/names input and a new tree is selected, fill 
    //the information in the new table to represent the information in the new tree.
    for(var i=0; i<this._ssm_ids.length; i++)
      this._add_row(num_samples,i);
    
    for(var i=0; i<this._ssm_ids.length; i++)
      this._fill_row(tidx,num_samples,i);
  }

  this._add_row = function(num_samples, row_idx){
    var new_row = this._first_row.clone().appendTo(this._table.find(".ssm-table-body"));
    new_row.find('.ssm-input').attr('style','display: none');
    new_row.find('#ssms-name-0').attr('id','ssms-name-' + (row_idx+1)).text('');
    new_row.find('#ssms-id-0').attr('id','ssms-id-' + (row_idx+1)).text('');
    new_row.find('#ssms-node-0').attr('id','ssms-node-' + (row_idx+1)).text('');
    for(var i=0;i<num_samples;i++)
      new_row.find('#vaf-samp' + i + '-0').attr('id','vaf-samp' + i + '-' + (row_idx+1)).text('');
  }

  this._fill_row = function(tidx, num_samples, row_idx){
    //With the information gathered, fill in the blank row in the table.
    this._input.val('');
    this._table.find('#ssms-id-' + row_idx).text(this._ssm_ids[row_idx]);
    this._table.find('#ssms-name-' + row_idx).text(this._ssm_names[row_idx]);
    //Get the read depth info from the muts file and fill in the vaf cells
    var tot_reads = this._muts.ssms[this._ssm_ids[row_idx]].total_reads;
    var ref_reads = this._muts.ssms[this._ssm_ids[row_idx]].ref_reads;
    for(var i=0;i<num_samples;i++){
      var vaf = (tot_reads[i] - ref_reads[i]) / tot_reads[i];
      var var_reads = tot_reads[i] - ref_reads[i];
      //entry text
      this._table.find('#vaf-samp' + i + '-' + row_idx).text(vaf.toFixed(3));
      //tooltip text
      this._table.find('#vaf-samp' + i + '-' + row_idx).attr("title","Variant Reads=" + var_reads + "; Total Reads=" + tot_reads[i]);
    }
    //Load in the mutass json, find which subpopulation the ssm belongs to, fill the cell
    var mutass_treepath = this._dataset.mutass_path + "/" + tidx  + '.json';
    var sViewer = this;
    d3.json(mutass_treepath, function(mutass){
      var pop = sViewer._find_pop_containing_ssm(mutass, sViewer._ssm_ids[row_idx]);
      sViewer._table.find('#ssms-node-' + row_idx).text(pop);
    })
  }
}

SSM_Viewer.prototype._find_pop_containing_ssm = function(mutass,ssm_id){
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

SSM_Viewer.prototype._set_ID_and_name_from_input_text = function(ssm_text){
  var success = false;
  if(ssm_text in this._muts.ssms){ //If typed in ID
    success = true;
    this._ssm_ids.push(ssm_text);
    if(this._muts.ssms[ssm_text].hasOwnProperty('name'))
      this._ssm_names.push(this._muts.ssms[ssm_text].name);
    else
      //insert a dash if the name isn't available
      this._ssm_names.push('\u2014');
  }else{ //If typed in nothing or name
    var sViewer = this;
    Object.keys(sViewer._muts.ssms).forEach(function(sidx){
      if(sViewer._muts.ssms[sidx].hasOwnProperty("name") & (sViewer._muts.ssms[sidx].name == ssm_text)){
        success = true;
        sViewer._ssm_ids.push(sidx);
        sViewer._ssm_names.push(ssm_text);
        return;
      }
    })
  }
  return success;
}