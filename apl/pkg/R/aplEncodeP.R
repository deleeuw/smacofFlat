#aplEncode <-
#function(ind,dims) {
#n<-length(dims); cell<-integer(n)
#if ((ind < 1) || (ind > prod(dims))) stop("No such cell")
#.C("aplEncodeC",as.integer(cell),as.integer(dims),as.integer(n),as.integer(ind))[[1]]
#}

aplEncode<-function(ind, dims){
    if( length( ind ) > 1 ){
       stop ( "Dimension error" )
    }
    if( ( ind < 1 ) || ( ind > prod( dims ) ) ){
       stop ( "No such cell" )
    }
    .Call( "APLENCODEP", as.integer( ind ), as.integer( dims ) )
}
