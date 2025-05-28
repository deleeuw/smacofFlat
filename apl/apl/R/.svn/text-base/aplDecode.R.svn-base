aplDecode <- function( cell, dims ){
    if( length( cell ) != length( dims ) ){
       stop( "Dimension error" )
    }
    if( any(cell > dims ) || any ( cell < 1 ) ){
       stop( "No such cell" )
    }
    .Call( "APLDECODE", as.integer( cell ), as.integer( dims ) )
}
