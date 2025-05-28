aplInnerProduct <- function( a, b, f = "*", g = "+" ){
    sa <- aplShape( a ); ra<- aplRank(a); ia <- 1:( ra - 1 ); 
    sb <- aplShape( b ); rb<- aplRank(b); ib <- ( ra - 1 ) + ( 1:( rb - 1 ) );
    if( !is.function( f ) ){ f <- match.fun( f ); }
    if( !is.function( g ) ){ g <- match.fun( g ); }
    env <- new.env( );
    environment( f ) <- env; environment( g ) <- env;
    ns <- last( sa ); nt <- first( sb );
    if ( ns != nt ) { stop( "Incompatible array dimensions" ) }
    sz <- c( butLast( sa ), butFirst( sb ) ); 
    z  <- .Call( "APLINNERPRODUCT",
                  f,
                  g,
                  as.double( a ),
                  as.double( b ),
                  as.integer( sa ),
                  as.integer( sb ),
                  as.integer( sz ),
                  as.integer( ns ),
                  env 
               );
    return(array( z, sz ) )
}

