#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SQUARE(x) ((x) * (x))

void smacofSSUREngine(int *nobj, int *ndim, int *ndat, int *itel, int *itmax,
                      int *digits, int *width, bool *verbose, double *sold,
                      double *snew, double *eps, int *iind, int *jind,
                      double *edis, double *dhat, double *xold, double *xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    while (true) {
        for (int k = 0; k < Ndat; k++) {
            int is = iind[k] - 1, js = jind[k] - 1;
            double elem = dhat[k] / edis[k];
            for (int s = 0; s < Ndim; s++) {
                double add = elem * (xold[is] - xold[js]);
                xnew[is] += add;
                xnew[js] -= add;
                is += Nobj;
                js += Nobj;
            }
        }
        for (int i = 0; i < Nobj; i++) {
            int is = i;
            for (int s = 0; s < Ndim; s++) {
                xnew[is] = xnew[is] / (double)Nobj;
                is += Nobj;
            }
        }
        for (int k = 0; k < Ndat; k++) {
            int is = iind[k] - 1, js = jind[k] - 1;
            double sum = 0.0;
            for (int s = 0; s < Ndim; s++) {
                sum += SQUARE(xnew[is] - xnew[js]);
                edis[k] = sqrt(sum);
                is += Nobj;
                js += Nobj;
            }
        }
        *snew = 0.0;
        for (int k = 0; k < Ndat; k++) {
            *snew += SQUARE(dhat[k] - edis[k]);
        }
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