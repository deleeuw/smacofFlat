aplRotate <- function( a, x, k ) 
{
    if( is.vector( a ) ) { return( aplRTV( a, k ) );}
    sa <- aplShape( a ); sx <- aplShape( x );
    if( aplRank( x ) == 1 ){ x <- array( x, sa[-k] ); }
    if( !identical( sa[-k], aplShape( x ) ) ){
        stop("Index Error")
    }
    z <- array( 0, sa ); sz <- sa; nz <- prod( sz ); sk<-sz[k]
    for( i in 1:nz ){
        ivec <- aplEncode( i, sz );
        xx <- x[aplDecode( ivec[-k], sx )];
        ak <- rep( 0, sk );
        for( j in 1:sk ){
            jvec <- ivec; jvec[k] <- j;
            ak[j] <- a[aplDecode( jvec, sz )]
            }       
        bk <- aplRTV( ak, xx )
        for( j in 1:sk ){
            jvec <- ivec; jvec[k] <- j;
            z[aplDecode( jvec, sz )] <- bk[j];
        }
    }
return( z );
}

aplRTV <- function( a, k ) 
{
    n <- aplShape( a )
    if( k > 0 )
        return( c( a[-(1:k)], a[1:k] ) )
    if( k < 0 )
        return( c( a[(n+k+1):n], a[1:(n+k)] ) )
    return( a )
}

