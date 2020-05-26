$(document).ready(function() {
    var table = $('#bigtable').DataTable( {
		"scrollX": true,
		"lengthMenu": [ [10, 25, 50, -1], [10, 25, 50, "All"] ],
    "stateSave": true,
    "stateDuration": 31556952,
    buttons: ['copy','colvis',
                {
                  text:'Hide Similar Columns',
                  action: function()
                  {
                    table.columns().every( function ( ) {
                      var data = this.data().unique();
                      if(data.length == 1)
                      {
                        console.log(data.length);
                        this.visible(false);
                      }
                    } );
                  }
                
                },
                {
                  text:'Show All Columns',
                  action: function()
                  {
                    table.columns().every( function ( ) {
                      this.visible(true);
                    } );
                  }
                
                }
              ],
    } );
    table.buttons().container()
        .appendTo( '#bigtable_wrapper .col-sm-6:eq(0)' );
});
