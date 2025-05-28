#include <math.h>
#include <stdio.h>

#define SQUARE(x) ((x) * (x))

void smacofDistances(int *nobj, int *ndim, int *ndat, int *ii, int *jj,
                     double *d, double *x) {
    int Nobj = *nobj, Ndim = *ndim, Ndat = *ndat;
    for (int k = 0; k < Ndat; k++) {
        int is = ii[k] - 1, js = jj[k] - 1;
        double sum = 0.0;
        for (int s = 0; s < Ndim; s++) {
            sum += SQUARE(x[is] - x[js]);
            d[k] = sqrt(sum);
            is += Nobj;
            js += Nobj;
        }
    }
}

void smacofStepUnweighted(int *nobj, int *ndim, int *ndat, int *ii, int *jj,
                          double *edis, double *dhat, double *xold,
                          double *xnew) {
    int Nobj = *nobj, Ndim = *ndim, Ndat = *ndat;
    for (int k = 0; k < Ndat; k++) {
        int is = ii[k] - 1, js = jj[k] - 1;
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
            printf("%f %f\n", xold[is], xnew[is]);
        }
    }
    for (int k = 0; k < Ndat; k++) {
      int is = ii[k] - 1, js = jj[k] - 1;
      double sum = 0.0;
      for (int s = 0; s < Ndim; s++) {
        sum += SQUARE(xnew[is] - xnew[js]);
        edis[k] = sqrt(sum);
        is += Nobj;
        js += Nobj;
      }
    }    
    return;
}

void smacofStepWeighted(int *nobj, int *ndim, int *ndat, int *ii, int *jj,
                        double *edis, double *dhat, double *wght, double *vinv,
                        double *xold, double *xtmp, double *xnew) {
    int Nobj = *nobj, Ndim = *ndim, Ndat = *ndat;
    for (int k = 0; k < Ndat; k++) {
        int is = ii[k] - 1, js = jj[k] - 1;
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
    (void)smacofDistances(nobj, ndim, ndat, ii, jj, edis, xnew);
    return;
}
