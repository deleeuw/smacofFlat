#include <R.h>
#include <Rinternals.h>
#include <math.h>
//#include <omp.h>

static inline void apldecode( int*, int*, int*, int* );
static inline void aplencode( int*, int*, int*, int*); 
SEXP APLDECODE( SEXP, SEXP );
SEXP APLENCODE( SEXP, SEXP );
SEXP APLSELECT( SEXP, SEXP, SEXP );
SEXP APLTRANSPOSE( SEXP, SEXP, SEXP, SEXP, SEXP );
SEXP APLSCAN( SEXP, SEXP, SEXP, SEXP, SEXP );
SEXP APLREDUCE( SEXP, SEXP, SEXP, SEXP, SEXP, SEXP );
SEXP APLINNERPRODUCT( SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP );

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

SEXP
APLDECODEP( SEXP cell, SEXP dims )
{
    int n = length( dims ), nProtected = 0; 
    int DIMS[n], CELL[n], IND;
    SEXP ind;
    #pragma omp parallel for
    for( int i = 0; i < n; i++ ){
       DIMS[i] = INTEGER( dims )[i];
       CELL[i] = INTEGER( cell )[i];
    }
    PROTECT( ind = allocVector( INTSXP, 1 ) ); ++nProtected;
    INTEGER( ind )[0] = 1;
    (void) apldecode (CELL, DIMS, &n, &IND);
    INTEGER(ind)[0] = IND;
    UNPROTECT( nProtected );
    return ind;
} 

SEXP
APLENCODEP( SEXP ind, SEXP dims )
{
    int  n = length( dims ), aux = INTEGER( ind )[0], nProtected = 0;
    int  CELL[n], DIMS[n];
    SEXP cell;
    PROTECT( cell = allocVector( INTSXP, n ) ); ++nProtected;
    #pragma omp parallel for
    for( int i = 0; i < n; i++ ){
       DIMS[i] = INTEGER( dims )[i];
    }
    (void) aplencode( CELL, DIMS, &n, &aux ); 
    #pragma omp parallel for
    for( int i = 0; i < n; i++ ){
       INTEGER(cell)[i] = CELL[i];
    }
    UNPROTECT( nProtected );
    return cell;
}
SEXP
APLSELECTP( SEXP a, SEXP dima, SEXP list )
{
   int  r = length( dima ), lz = 1, dimzi, nind, itel, czll[r], cell[r], dimz[r], da[r], nProtected = 0;
   SEXP   z;

   #pragma omp parallel for reduction(*:lz) private(dimzi)
   for( int i = 0; i < r; i++ ){
       dimzi = length( VECTOR_ELT( list, i ) );
       dimz[i] = dimzi;
       lz *= dimzi;
       da[i] = INTEGER( dima )[i];
   }
   PROTECT( z = allocVector( REALSXP, lz ) );  ++nProtected;
   for( int i = 0; i < lz; i++ ){
       itel = i + 1;
       (void) aplencode( cell, dimz, &r, &itel ); 
       for ( int j = 0; j < r; j++ ) {
           czll[j] = INTEGER ( VECTOR_ELT( list, j ) )[cell[j] - 1];
       }
       (void) apldecode (czll, da, &r, &nind);
       REAL( z )[i] = REAL( a )[nind - 1];
   }
   UNPROTECT( nProtected );
   return( z );
}

SEXP
APLTRANSPOSEP( SEXP a, SEXP x, SEXP sa, SEXP sz, SEXP rz )
{
    int  itel, nind, na = 1, nz = 1, ra = length( sa ), lsz = length( sz ), nProtected=0;
    int  RZ = INTEGER(rz)[0];
    int  ivec[RZ], jvec[ra], SA[ra], SZ[lsz];
    SEXP z;
    #pragma omp parallel for reduction(*:na)
    for( int i = 0; i < ra ; i++ ){
        SA[i] = INTEGER(sa)[i];
        na *= INTEGER( sa )[i]; 
    }
    #pragma omp parallel for reduction(*:nz)
    for( int i = 0; i < lsz; i++ ){ 
        SZ[i] = INTEGER( sz )[i]; 
        nz *= INTEGER( sz )[i]; 
    }
    PROTECT( z    = allocVector( REALSXP,            nz  ) ); ++nProtected;
    for( int i = 0; i < nz; i++ ){
        itel = i + 1;
        (void) aplencode( ivec, SZ, &RZ, &itel);
        for( int j = 0; j < ra; j++ ){
            jvec[j] = ivec[INTEGER( x )[j] - 1];
        }
        (void) apldecode(jvec, SA, &ra, &nind);
        REAL( z )[i] = REAL( a )[nind - 1];
    }
    UNPROTECT( nProtected );
    return z;
}

SEXP
APLREDUCEP( SEXP f, SEXP a, SEXP k, SEXP sa, SEXP sz, SEXP env )
{
    int u, r, kk, na = 1, nz = 1, nProtected = 0;
    int nk = length( k ), ra = length( sa ), rz = length( sz );
    int SA[ra], SZ[rz];
    SEXP z, R_fcall= R_NilValue;
    SEXP Z = R_NilValue, A = R_NilValue;
    #pragma omp parallel for reduction(*:na)
    for( int i = 0; i < ra; i++ ){ 
        SA[i] = INTEGER( sa )[i]; 
        na *= INTEGER( sa )[i]; 
    }
    #pragma omp parallel for reduction(*:nz)
    for( int i = 0; i < rz; i++ ){ 
        SZ[i] = INTEGER( sz )[i]; 
        nz *= INTEGER( sz )[i]; 
    }
    kk = ra - nk;
    int itel, nind, ivec[ra],kvec[kk], ind[nz];
    PROTECT( R_fcall= lang3(f, R_NilValue, R_NilValue) ); ++nProtected;
    PROTECT( Z      = allocVector( REALSXP,       1  ) ); ++nProtected;
    PROTECT( A      = allocVector( REALSXP,       1  ) ); ++nProtected;
    PROTECT( z      = allocVector( REALSXP,       nz ) ); ++nProtected;
    #pragma omp parallel for
    for( int i = 0; i < nz; i++ ){
        ind[i] = 0;
    }
    for( int i = 0; i < na; i++ ){
        itel = i + 1;
        (void) aplencode(ivec,SA,&ra,&itel);
        u = 0;
        for( int j = 0; j < ra; j++ ){
            r = 0;
            for( int v = 0; v < nk; v++ ){
                if( j == ( INTEGER( k )[v] - 1 ) ) r = 1;
            }
            if( r == 0 ){
                kvec[u] = ivec[j];
                u += 1;
            }
        }
        (void) apldecode(kvec,SZ,&rz,&nind);
        if ( ind[nind - 1] == 0 ) {
            REAL( z )[nind - 1] = REAL( a )[i];
            ind[nind - 1] = 1;
        } else {
            REAL( Z )[0] = REAL( z )[nind - 1];
            REAL( A )[0] = REAL( a )[i];
            SETCADR( R_fcall, Z );
            SETCADDR( R_fcall, A );
            REAL( z )[nind - 1] = REAL( eval( R_fcall, env ) )[0];
        }
    }
    UNPROTECT( nProtected );
    return z;
}

SEXP 
APLSCANP( SEXP f, SEXP a, SEXP k, SEXP sa, SEXP env )
{
    int sk, l, na=1, itel, nind, ra = length( sa ),nProtected=0;
    int ivec[ra], SA[ra];
    #pragma omp parallel for reduction(*:na)
    for( int i = 0; i < ra; i++ ){ 
        SA[i]=INTEGER( sa )[i];
        na *= INTEGER( sa )[i]; 
    }
    SEXP z,  Z, A, R_fcall = R_NilValue;
    PROTECT( R_fcall = lang3( f, R_NilValue, R_NilValue ) ); ++nProtected;
    PROTECT( Z       = allocVector( REALSXP,         1  ) ); ++nProtected;
    PROTECT( A       = allocVector( REALSXP,         1  ) ); ++nProtected;
    PROTECT( z       = allocVector( REALSXP,         na ) ); ++nProtected;
    l = INTEGER( k )[0] - 1;
    for( int i = 0; i < na; i++ ){
        itel = i + 1;
        (void) aplencode( ivec, SA, &ra, &itel );
        sk = ivec[l];
        if( sk == 1 ){
             REAL( z )[i] = REAL( a )[i];
        }else{
            ivec[l] -= 1;
            (void) apldecode( ivec, SA, &ra, &nind );
            REAL( Z )[0]=REAL( z )[nind-1];
            REAL( A )[0]=REAL( a )[i];
            SETCADR( R_fcall, Z );
            SETCADDR( R_fcall, A );
            REAL( z )[i] = REAL( eval( R_fcall, env ) )[0];
        }
    }
    UNPROTECT( nProtected );
    return z;
}



SEXP
APLINNERPRODUCTP(SEXP f, SEXP g, SEXP a, SEXP b, SEXP sa, SEXP sb, SEXP sz, SEXP ns, SEXP env)
{
    int nz = 1, nProtected = 0;
    int ra = length( sa ),rb = length( sb ), rz = length( sz );
    int SZ[rz], SA[ra], SB[rb];
    SEXP z, A, B, Z, t;
    SEXP R_fcall = R_NilValue, R_gcall = R_NilValue;
    #pragma omp parallel for reduction(*:nz)
    for( int i = 0; i < rz; i++ ){ 
        SZ[i] = INTEGER( sz )[i]; 
        nz *= INTEGER( sz )[i]; 
    }
    #pragma omp parallel for
    for( int i = 0; i < ra; i++ ){ 
        SA[i] = INTEGER( sa )[i]; 
    }
    #pragma omp parallel for
    for( int i = 0; i < rb; i++ ){ 
        SB[i] = INTEGER( sb )[i]; 
    }
    int ivec[rz], jvec[ra], kvec[rb], itel, k, l;
    PROTECT( R_fcall = lang3( f, R_NilValue, R_NilValue ) ); ++nProtected;
    PROTECT( R_gcall = lang3( g, R_NilValue, R_NilValue ) ); ++nProtected;
    PROTECT( z       = allocVector( REALSXP,         nz ) ); ++nProtected;
    PROTECT( t       = allocVector( REALSXP,         1  ) ); ++nProtected;
    PROTECT( Z       = allocVector( REALSXP,         1  ) ); ++nProtected;
    PROTECT( B       = allocVector( REALSXP,         1  ) ); ++nProtected;
    PROTECT( A       = allocVector( REALSXP,         1  ) ); ++nProtected;
    for( int i = 0; i < nz; i++ ){
        itel = i + 1;
        (void) aplencode(ivec, SZ, &rz,&itel);
        for( int j = 0; j < INTEGER( ns )[0]; j++ ){
            #pragma omp parallel for
            for( int u = 0; u < ra - 1; u++ ){
                jvec[u] = ivec[u];
            }
            jvec[ra - 1] = j + 1;
            (void) apldecode( jvec, SA, &ra, &k);
            for( int u = 1; u < rb; u++ ){
                kvec[u] = ivec[ra + u - 2];
            }
            kvec[0] = j + 1;
            (void) apldecode(kvec, SB, &rb, &l);
            REAL( A )[0] = REAL( a )[k - 1];
            REAL( B )[0] = REAL( b )[l - 1];
            SETCADR( R_fcall, A );
            SETCADDR( R_fcall, B );
            REAL( t )[0] = REAL( eval( R_fcall, env ) )[0];
            if( j == 0 ){ 
                REAL( z )[i] = REAL( t )[0];
            } else {
                REAL( Z )[0] = REAL( z )[i];
                SETCADR( R_gcall, t );
                SETCADDR( R_gcall, Z );
                REAL( z )[i] = REAL( eval( R_gcall, env ) )[0];
            }
        }
    }
    UNPROTECT( nProtected );
    return z;
}

