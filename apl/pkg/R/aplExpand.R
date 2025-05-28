aplExpand <- function( x, y, axis = 1 ){
    if( is.vector( x ) ){ return( aplEXV( x, y ) ) }
    d <- dim( x ); m <- which( y ); e <- d; e[axis] <- m;
    if( m != d[axis] ){
        stop("Incorrect dimension length")
    }
    z <- array( 0, e )
    apply( z,( 1:n )[-axis],function(i) z[i]<-x[i])
}

aplEXV <- function( x, y ){
    z <- rep( 0, length( y ) )
    m <- which( y )
    if( length( m ) != length( x ) ){
        stop("Incorrect vector length")
    }
    z[which( y )] <- x
    return( z )
}

