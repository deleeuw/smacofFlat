aplGet <- function( a, cell ) 
{
    dims <- dim( a ); n <- length( dims ); b <- 0;
    if( any( cell > dims ) || any( cell < 1 ) ){ stop("No such cell"); }
    return( a[aplDecode( cell, dims )]);
} 

