aplReduceP <- function( a, k, f = "+" ){
    if( is.vector( a ) ){
        return( aplRDV( a, f ) );
    }
    nk  <- length(k);
    sa  <- aplShape(a); 
    sz  <- sa[(1:length(sa))[-k]]; 
    if( !is.function( f ) ){ f <- match.fun( f ); }
    env <- new.env( );
    environment( f ) <- env;
    z <- .Call( "APLREDUCEP", 
                f, 
                as.double(a), 
                as.integer(k),  
                as.integer(sa), 
                as.integer(sz), 
                env
              );
    return( array( z, sz ) )
}

aplRDV <- function( x, f = "+" ) 
{
    if ( length( x ) == 0 ){ return( x ); }
    s <- x[1]
    if( length( x ) == 1 ){ return(s) }
    for ( i in 2:length( x ) ){
        s <- match.fun( f )( s, x[i] )
    }
    return(s)
}
