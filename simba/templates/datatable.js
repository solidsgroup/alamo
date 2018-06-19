$(document).ready(function() {
    var table = $('#bigtable').DataTable( {
		"scrollX": true,
		"lengthMenu": [ [10, 25, 50, -1], [10, 25, 50, "All"] ],
		stateSave: true,
		buttons: ['copy','colvis'],
    } );
    table.buttons().container()
        .appendTo( '#bigtable_wrapper .col-sm-6:eq(0)' );
});


