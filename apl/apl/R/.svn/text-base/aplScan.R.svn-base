aplScan <- function( a, k, f = "+" ){
    if( is.vector( a ) ){ return( aplSCV( a, f ) ); }
    if( !is.function( f ) ){ f <- match.fun( f ); }
    env <- new.env( );
    environment( f ) <- env;
    sa <-aplShape(a); ra <-aplRank(a); 
    z  <- .Call( "APLSCAN",
                 f, 
                 as.double( a ), 
                 as.integer( k ),
                 as.integer( sa ),
                 as.integer( ra ),
                 env
               );
    return( array( z, sa ) );
}

aplSCV <- function( x, f = "+" ){
    if ( length( x ) <= 1 ) return( x )
    return( sapply( 1:length( x ), function( i ) aplRDV( x[1:i], f ) ) )
}
