#include <math.h>
#include <stdio.h>

#define SQUARE(x) ((x) * (x))

void smacofDistances(int *nobj, int *ndim, int *ndat, int *ii, int *jj,
                     double *edis, double *xnew) {
    int Nobj = *nobj, Ndim = *ndim, Ndat = *ndat;
    for (int k = 0; k < Ndat; k++) {
        int i = ii[k], j = jj[k];
        double sum = 0.0;
        for (int s = 0; s < Ndim; s++) {
            int is = s * Nobj + (i - 1), js = s * Nobj + (j - 1);
            sum += SQUARE(xnew[is] - xnew[js]);
            edis[k] = sqrt(sum);
        }
    }
}

void smacofStepUnweighted(int *nobj, int *ndim, int *ndat, int *ii, int *jj,
                          double *edis, double *dhat, double *xold,
                          double *xnew) {
    int Nobj = *nobj, Ndim = *ndim, Ndat = *ndat;
    for (int k = 0; k < Ndat; k++) {
        int i = ii[k], j = jj[k];
        double elem = dhat[k] / edis[k];
        for (int s = 0; s < Ndim; s++) {
            int is = s * Nobj + (i - 1), js = s * Nobj + (j - 1);
            double add = elem * (xold[is] - xold[js]);
            xnew[is] += add;
            xnew[js] -= add;
        }
    }
    for (int i = 0; i < Nobj; i++) {
        for (int s = 0; s < Ndim; s++) {
            int is = s * Nobj + i;
            xnew[is] = xnew[is] / (double)Nobj;
        }
    }
    (void)smacofDistances(nobj, ndim, ndat, ii, jj, edis, xnew);
    return;
}

void smacofStepWeighted(int *nobj, int *ndim, int *ndat, int *ii, int *jj,
                        double *edis, double *dhat, double *wght, double *vinv,
                        double *xold, double *xtmp, double *xnew) {
    int Nobj = *nobj, Ndim = *ndim, Ndat = *ndat;
    for (int k = 0; k < Ndat; k++) {
        int i = ii[k], j = jj[k];
        double elem = wght[k] * dhat[k] / edis[k];
        for (int s = 0; s < Ndim; s++) {
            int is = s * Nobj + (i - 1), js = s * Nobj + (j - 1);
            double add = elem * (xold[is] - xold[js]);
            xtmp[is] += add;
            xtmp[js] -= add;
        }
    }
    int k = 0;
    for (int j = 0; j < Nobj - 1; j++) {
        for (int i = j + 1; i < Nobj; i++) {
            double elem = vinv[k];
            for (int s = 0; s < Ndim; s++) {
                int is = s * Nobj + i, js = s * Nobj + j;
                double add = elem * (xtmp[is] - xtmp[js]);
                xnew[is] += add;
                xnew[js] -= add;
            }
            k++;
        }
    }
    (void)smacofDistances(nobj, ndim, ndat, ii, jj, edis, xnew);
    return;
}
