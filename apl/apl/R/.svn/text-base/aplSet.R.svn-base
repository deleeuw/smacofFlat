aplSet <-
function(a,b,cell) {
dims<-dim(a); n<-length(dims)
if (any(cell>dims) || any(cell<1)) stop("No such cell")
a[aplDecode(cell,dims)]<-b
return(a)
}

