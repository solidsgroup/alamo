$(document).ready(function() {
    var table = $('#bigtable').DataTable( {
		"scrollX":true
    } );
	
    $('input.toggle-vis').on( 'change', function (e) {
        e.preventDefault();
        var column = table.column( $(this).attr('data-column') );
        //column.visible( ! column.visible() );
		column.visible( $(this).is(":checked") );
		table.columns.adjust().draw();
		$('input.toggle-vis').addClass('btn-xs');
    } );


  $('.all').on('click', function(e){
    $this = this;  
    $.each($(this).parents('ul').find('input'), function(i, item){
		$(item).prop('checked', $this.checked);
		$(item).trigger('change');
    });
  });

} );

