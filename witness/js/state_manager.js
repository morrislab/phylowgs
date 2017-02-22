var StateManager = {
  _state: {},

  // Count number of restore events currently in progress.
  restoring: 0,

  update: function(key, value) {
    if(this.restoring > 0)
      return;
    if(value == null) {
      delete this._state[key];
    } else {
      this._state[key] = value;
    }

    var serialized = btoa(JSON.stringify(this._state));
    window.location.hash = serialized;
  },

  _fetch_state: function() {
    var state = window.location.hash.substring(1);
    if(state.length == 0)
      return null;
    state = JSON.parse(atob(state));
    return state;
  },

  restore: function() {
    var self = this;
    var state = this._fetch_state();
    if(state === null)
      return;
    this._state = state;
    // Don't let state elements be overwritten when we call the click() method
    // of various elements.

    if(state.hasOwnProperty('section')) {
      self.restoring++;
      $('.navbar-nav a').filter('.' + state.section).click();
      self.restoring--;
    }

    if(state.hasOwnProperty('run')) {
      self.restoring++;
      $('#runs a').filter(function() {
        return $(this).text() === state.run;
      }).click();
      self.restoring--;
    }

    if(state.hasOwnProperty('sample')) {
      self.restoring++;
      $('#samples').find('a').filter(function() {
        return $(this).text() === state.sample;
      }).click();
      self.restoring--;
    }

    // I don't know why I need to delay this but not the other state-restoring
    // parts. Meh. Whatever. This is all a giant hack, anyway.
    if(state.hasOwnProperty('tidx')) {
      self.restoring++;
      setTimeout(function() {
        $('#trees tbody tr').filter(function() {
          return $(this).find('.tree-index').text() === state.tidx;
        }).click();
        self.restoring--;
      }, 500);
    }
  }
};
