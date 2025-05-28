aplTake <- function(a,x,drop=FALSE) {
    sa<- aplShape( a ); ra<-aplRank( a );
    y <- as.list( rep( 0, ra ) );
    if( ra != length( x ) ){
       stop( "Dimension error" );
    }
    if( any(sa <= abs( x ) ) ){
       stop( "Cannot drop more than the maximum dimension." );
    }
    if( any( x == 0 ) ){
        return( NULL);
    }
    for (i in 1:ra) {
        ss <- sa[i]; xx <- x[i]; sx <- ss-xx;
        if ( xx > 0 ) y[[i]]<-1:xx
        if ( xx < 0 ) y[[i]]<-(sx+1):ss
        }
    return( aplSelect( a, y, drop ) );
}

