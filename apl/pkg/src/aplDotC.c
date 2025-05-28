#include <R.h>
#include <Rinternals.h>

static char* f2_char;
static char* g2_char;

void aplDecodeC(int* cell, int* dims, int* n, int* ind) {
int i, aux = 1; *ind = 1;
for (i = 0; i < *n; i++) {
    *ind += aux * (cell[i] - 1);
    aux *= dims[i];
    }
}

void aplEncodeC(int* cell, int* dims, int* n, int* ind) {
int i, aux = *ind, nn = *n, pdim = 1;
for (i = 0; i < nn - 1; i++)
    pdim *= dims[i];
for (i = nn - 1; i > 0; i--) {
    cell[i] = (aux - 1)/pdim;
    aux -= pdim*cell[i];
    pdim /= dims[i-1];
    cell[i] += 1;
    }
cell[0] = aux;
}

void aplSelectC(double *a, int *sa, int *ra, SEXP *list, double *z, int *sz, int *rz, int *nz)
{
int i, j, k; 
int ivec[*rz], jvec[*ra];
for (i = 0; i < *nz; i++) {
	k = i + 1;
	(void) aplEncodeC (ivec,sz,rz,&k); 
    for (j = 0; j < *ra; j++) {
        SEXP vec = list[j];
        jvec[j] = INTEGER(vec)[ivec[j]-1];
    }
	k = 1;
	(void) aplDecodeC (jvec,sa,ra,&k);
	z[i] = a[k - 1];
    }
}

void aplTransposeC(double *a, int *x, int *sa, int *ra, int *na, int *sz, int *rz, int *nz, double *z)
{
int i, j, r, ivec[*rz], jvec[*ra];
for (i = 0; i < *nz; i++){
	r = i + 1;
	(void) aplEncodeC(ivec,sz,rz,&r);
	for (j = 0; j < *ra; j++)
		jvec[j] = ivec[x[j]-1];
	(void) aplDecodeC(jvec,sa,ra,&r);
	z[i] = a[r - 1];
	}
}

void aplReduceC(char** funclist, double *a, int *k, int *nk, int *na, int *sa, int *ra, int *nz, int *sz, int *rz, double *z) {
double f2_glue(double, double);
double f2_comp(double (*)(),double,double);
int i, j, u, v, r, m, kk = (*ra) - (*nk), ivec[*ra], kvec[kk], ind[*nz];
for (i = 0; i < *nz; i++) ind[i] = 0;
f2_char = funclist[0];
for (i = 0; i < *na; i++){
    r = i + 1;
	(void) aplEncodeC(ivec,sa,ra,&r);
	u = 0;
	for (j = 0; j < *ra; j++) {
		r = 0;
		for (v = 0; v < *nk; v++) {
			if (j == (k[v] - 1)) r = 1;
		    }
		if (r == 0)	{
			kvec[u] = ivec[j];
			u += 1;
			}
		}
	(void) aplDecodeC(kvec,sz,rz,&m);
	if (ind[m - 1] == 0) {
		z[m - 1] = a[i];
		ind[m - 1] = 1;
		}
	else
		z[m - 1] = f2_comp((double(*)())f2_glue,z[m - 1],a[i]);
	}
}

void aplScanC(char** funclist,double *a, int *k, int *na, int *sa, int *ra, double *z)
{
double f2_glue(double, double);
double f2_comp(double (*)(),double,double);
int i, r, sk, l, ivec[*ra];
l = *k - 1;
f2_char = funclist[0];	
for (i = 0; i < *na; i++){
    r = i + 1;
	(void) aplEncodeC(ivec,sa,ra,&r);
	sk = ivec[l];
	if (sk == 1) z[i] = a[i];
		else {
			ivec[l] -= 1;
			(void) aplDecodeC(ivec,sa,ra,&r);
			z[i] = f2_comp((double(*)())f2_glue,z[r - 1],a[i]);
		    }
    }
}

void aplInnerProductC(char** funclist, double *a, double *b, int *sa, int *ra, int *sb, int *rb, int *sz, int *rz, int *nz, int *ns, double *z) {
double f2_glue(double, double);
double g2_glue(double, double);
double f2_comp(double (*)(),double,double);
double g2_comp(double (*)(),double,double);
double t; int i, j, r, k, l, u, ivec[*rz], jvec[*ra], kvec[*rb];
f2_char = funclist[0]; g2_char= funclist[1]; k = l = 0;
for (i = 0; i < *nz; i++) {
	r = i + 1;
	(void) aplEncodeC(ivec,sz,rz,&r);
	for (j = 0; j < *ns; j++) {
		for (u = 0; u < *ra - 1; u++)
			jvec[u] = ivec[u];
		jvec[*ra - 1] = j + 1;
		(void) aplDecodeC(jvec,sa,ra,&k);
		for (u = 1; u < *rb; u++)
			kvec[u] = ivec[*ra + u - 2];
		kvec[0] = j + 1;
		(void) aplDecodeC(kvec,sb,rb,&l);
		t = f2_comp((double(*)())f2_glue,a[k-1],b[l-1]);
		if (j == 0) z[i] = t;
	        else z[i] = g2_comp((double(*)())g2_glue,t,z[i]);
    	}
	}
}

double f2_glue(double x, double y){ 
char *modes[2], *values[1];
void *args[2]; 
double xx[1], yy[1], *result; 
long lengths[2], nargs = (long)2,  nvals = (long)1;
lengths[0] = lengths[1] = (long)1;
nargs = (long)2;					 
args[0] = (void *)xx; xx[0] = x;
args[1] = (void *)yy; yy[0] = y;					 
modes[0] = modes[1] = "double"; 
call_R(f2_char,nargs,args,modes,lengths,0,nvals,values); 
result = (double*)values[0]; 
return(result[0]); 
} 

double g2_glue(double x, double y){ 
char *modes[2], *values[1];
void *args[2]; 
double xx[1], yy[1], *result; 
long lengths[2], nargs = (long)2,  nvals = (long)1;
lengths[0] = lengths[1] = (long)1;
nargs = (long)2;					 
args[0] = (void *)xx; xx[0] = x;
args[1] = (void *)yy; yy[0] = y;					 
modes[0] = modes[1] = "double"; 
call_R(g2_char,nargs,args,modes,lengths,0,nvals,values); 
result = (double*)values[0]; 
return(result[0]); 
} 

double f2_comp(double (*f)(), double x, double y) {
return f(x,y);
}

double g2_comp(double (*g)(), double x, double y) {
return g(x,y);
}

