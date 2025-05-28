aplSelectP <- function ( a, x, drop = TRUE ){
    dima = aplShape( a );
    if( length( dima ) != length( x ) ) {
        stop( "Dimension error" );
    }
    z <- .Call( "APLSELECTP", 
                as.double( a ), 
                as.integer( dima ), 
                lapply( x, as.integer ) 
              );
    z <- array( z, sapply( x, length ) );
    if( drop ){
        return( drop( z ) );
    }
    return( z );
}
