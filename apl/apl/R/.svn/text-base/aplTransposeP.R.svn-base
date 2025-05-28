aplTransposeP <- function( a, x = rev( 1:aplRank( a ) ) ) {
    sa <- aplShape( a ); 
    if( length( x ) != length( sa ) ) { stop( "Length Error" ) }
    rz <- max( x ); sz <- rep( 0, rz )
    for( i in 1:rz ){
        sz[i] <- min( sa[which( x == i )] )
    }
    z <- .Call( "APLTRANSPOSEP",
                as.double(a),
                as.integer( x ),
                as.integer( sa ),
                as.integer( sz ),
                as.integer( rz )
              );
    return( array( z, sz ) );
}

