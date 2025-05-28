aplDrop <- function( a, x, drop = FALSE )
{
    sa <- aplShape( a ); ra <- aplRank( a );
    y  <- as.list( rep( 0, ra ) );
    if( ra != length( x ) ){
       stop( "Dimension error" );
    }
    if( any(sa <= abs( x ) ) ){
       stop( "Cannot drop more than the maximum dimension." );
    }
    for( i in 1:ra ){
        ss <- sa[i]; 
        xx <- x[i]; 
        sx <- ss-xx
        if( xx >= 0 ){ y[[i]] <- ( xx + 1 ):ss; }
        if( xx <  0 ){ y[[i]] <- 1:sx; }
    }
    return( aplSelect( a, y, drop ) );
}
