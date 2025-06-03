#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SQUARE(x) ((x) * (x))

void smacofSSWREngine(int *nobj, int *ndim, int *ndat, int *itel, int *itmax,
                      int *digits, int *width, bool *verbose, double *wsum, double *sold,
                      double *snew, double *eps, int *iind, int *jind,
                      double *edis, double *dhat, double *wght, double *vinv,
                      double *xold, double *xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    while (true) {
        double *xtmp = (double *)calloc(Nobj * Ndim, sizeof(double));
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
        free(xtmp);
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
        *snew = 0.0;
        for (int k = 0; k < Ndat; k++) {
            *snew += wght[k] * SQUARE(dhat[k] - edis[k]);
        }
        *snew /= *wsum;
        if (*verbose) {
            printf("itel %4d sold %*.*f snew %*.*f\n", *itel, *width, *digits,
                   *sold, *width, *digits, *snew);
        }
        if ((*itel == *itmax) || ((*sold - *snew) < *eps)) {
            break;
        }
        for (int k = 0; k < Nobj * Ndim; k++) {
            xold[k] = xnew[k];
            xnew[k] = 0.0;
        }
        *sold = *snew;
        *itel += 1;
    }
}