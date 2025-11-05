#include "smacofSS.h"

void smacofSSUOEngine(int* nobj, int* ndim, int* ndat, int* itel, int* ties,
                      int* itmax, int* digits, int* width, bool* verbose,
                      double* sold, double* snew, double* eps, int* iind,
                      int* jind, int* blks, double* edis, double* dhat,
                      double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    while (true) {
        for (int k = 0; k < Nobj * Ndim; k++) {
            xnew[k] = 0.0;
        }
        for (int k = 0; k < Ndat; k++) {
            if (edis[k] == 0.0) {
                continue;
            }
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
        double smid = 0.0;
        for (int k = 0; k < Ndat; k++) {
            smid += SQUARE(dhat[k] - edis[k]);
        }
        smid /= (double)Ndat;
        double* wght = xcalloc(Ndat, sizeof(double));
        for (int k = 0; k < Ndat; k++) {
            wght[k] = 1.0;
        }
        dhat = memcpy(dhat, edis, Ndat * sizeof(double));
        if (*ties == 1) {
            (void)primaryApproach(ndat, blks, dhat, wght, edis, iind, jind);
        }
        if (*ties == 2) {
            (void)secondaryApproach(ndat, blks, dhat, wght);
        }
        if (*ties == 3) {
            (void)tertiaryApproach(ndat, blks, dhat, wght);
        }
        double ssq = 0.0;
        for (int k = 0; k < Ndat; k++) {
            ssq += SQUARE(dhat[k]);
        }
        *snew = 0.0;
        for (int k = 0; k < Ndat; k++) {
            dhat[k] *= sqrt(((double)Ndat) / ssq);
            *snew += SQUARE(dhat[k] - edis[k]);
        }
        *snew /= (double)Ndat;
        if (*verbose) {
            printf("itel %4d sold %*.*f smid %*.*f snew %*.*f\n", *itel, *width,
                   *digits, *sold, *width, *digits, smid, *width, *digits,
                   *snew);
        }
        if ((*itel == *itmax) || ((*sold - *snew) < *eps)) {
            break;
        }
        xold = memcpy(xold, xnew, Nobj * Ndim * sizeof(double));
        *sold = *snew;
        *itel += 1;
        xfree(wght);
    }
}
