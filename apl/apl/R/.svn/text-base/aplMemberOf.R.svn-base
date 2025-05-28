aplMemberOf <-
function(a,b) {
    if (!identical(typeof(a),typeof(b)))
        warning("Arguments of different types")
    arrTest(a); arrTest(b)
    sa<-aplShape(a); sb<-aplShape(b)
    na<-prod(sa); nb<-prod(sb)
    z<-array(0,sa)
    for (i in 1:na) {
        z[i]<-0; aa<-a[i]
        for (j in 1:nb)
            if (identical(aa,b[j])) z[i]<-1
        }
return(z)   
}

arrTest <-
function(x) if (!is.array(x)) stop("Functions in the apl package take array arguments")

