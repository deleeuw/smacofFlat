void aplReduceC(char** funclist, double *a, int *k, int *nk, int *na, int *sa, int *ra,
double f2_glue(double, double);
double f2_comp(double (*)(),double,double);
int i, j, u, v, r, m, kk = (*ra) - (*nk), ivec[*ra], kvec[kk], ind[*nz];
for (i = 0; i < *nz; i++) ind[i] = 0;
f2_char = funclist[0];
for (i = 0; i < *na; i++){
    r = i + 1;
    (void)aplEncodeC(ivec, sa, ra, &r);
    u = 0;
    for (j = 0; j < *ra; j++) {
        r = 0;
        for (v = 0; v < *nk; v++) {
            if (j == (k[v] - 1)) r = 1;
        }
        if (r == 0) {
            kvec[u] = ivec[j];
            u += 1;
        }
    }
    (void)aplDecodeC(kvec, sz, rz, &m);
    if (ind[m - 1] == 0) {
        z[m - 1] = a[i];
        ind[m - 1] = 1;
    } else
        z[m - 1] = f2_comp((double (*)())f2_glue, z[m - 1], a[i]);
}
46
}
{
void aplScanC(char**