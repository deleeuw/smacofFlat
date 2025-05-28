#include <R.h>
#include <Rinternals.h>
#include <math.h>

//static inline void apldecode( int*, int*, int*, int* );
//static inline void aplencode( int*, int*, int*, int*); 
static inline void
apldecode( int* cell, int* dims, int* n, int* ind )
{
    int  i, aux = 1; *ind = 1;
    for( i = 0; i < *n; i++ ) {
        *ind += aux * ( cell[i] - 1 );
        aux *= dims[i];
    }
} 

static inline void
aplencode( int* cell, int* dims, int* n, int* ind)
{
    int aux = *ind , pdim = 1;
    for( int i = 0; i < *n - 1; i++ ){
        pdim *= dims[i];
    }
    for( int i = *n - 1; i > 0; i-- ){
        cell[i] = ( aux - 1 ) / pdim;
        aux -= pdim * cell[i];
        pdim /= dims[i - 1];
        cell[i] += 1;
    }
    cell[0] = aux;
}

SEXP APLDECODE( SEXP, SEXP );
SEXP APLENCODE( SEXP, SEXP );
SEXP APLSELECT( SEXP, SEXP, SEXP );
SEXP APLTRANSPOSE( SEXP, SEXP, SEXP, SEXP, SEXP );
SEXP APLSCAN( SEXP, SEXP, SEXP, SEXP, SEXP );
SEXP APLREDUCE( SEXP, SEXP, SEXP, SEXP, SEXP, SEXP );
SEXP APLINNERPRODUCT( SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP );
