#include "smacofSS.h"

void smacofSSWMajorize(int* nobj, int* ndim, int* ndat, double* snew, int* iind,
                       int* jind, double* wght, double* vinv, double* edis,
                       double* dhat, double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    double* xtmp = xmalloc(Nobj * Ndim * sizeof(double));
    double wsum = 0.0;
    for (int k = 0; k < Nobj * Ndim; k++) {
        xtmp[k] = 0.0;
        wsum += wght[k];
    }
    for (int k = 0; k < Ndat; k++) {
        if (edis[k] == 0.0) {
            continue;
        }
        int is = iind[k] - 1, js = jind[k] - 1;
        double elem = wght[k] * dhat[k] / edis[k];
        for (int s = 0; s < Ndim; s++) {
            double add = elem * (xold[is] - xold[js]);
            xtmp[is] += add;
            xtmp[js] -= add;
            is += Nobj;
            js += Nobj;
        }
    }
    for (int k = 0; k < Nobj * Ndim; k++) {
        xnew[k] = 0.0;
    }
    int k = 0;
    for (int j = 0; j < Nobj - 1; j++) {
        for (int i = j + 1; i < Nobj; i++) {
            double elem = vinv[k];
            int is = i, js = j;
            for (int s = 0; s < Ndim; s++) {
                double add = elem * (xtmp[is] - xtmp[js]);
                xnew[is] += add;
                xnew[js] -= add;
                is += Nobj;
                js += Nobj;
            }
            k++;
        }
    }
    xfree(xtmp);
    for (int k = 0; k < Ndat; k++) {
        int is = iind[k] - 1, js = jind[k] - 1;
        double sum = 0.0;
        for (int s = 0; s < Ndim; s++) {
            sum += SQUARE(xnew[is] - xnew[js]);
            is += Nobj;
            js += Nobj;
        }
        edis[k] = sqrt(sum);
    }
    double sum = 0.0;
    for (int k = 0; k < Ndat; k++) {
        sum += wght[k] * SQUARE(dhat[k] - edis[k]);
    }
    *snew = sum / wsum;
    return;
}