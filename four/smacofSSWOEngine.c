#include "smacofSS.h"

void smacofSSWOEngine(int* nobj, int* ndim, int* ndat, int* itel, int* ties,
                      int* itmax, int* digits, int* width, bool* verbose,
                      double* wsum, double* sold, double* snew, double* eps,
                      int* iind, int* jind, int* blks, double* edis,
                      double* dhat, double* wght, double* vinv, double* xold,
                      double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    while (true) {
        double* xtmp = xmalloc(Nobj * Ndim * sizeof(double));
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
        double smid = 0.0;
        for (int k = 0; k < Ndat; k++) {
            smid += wght[k] * SQUARE(dhat[k] - edis[k]);
        }
        smid /= *wsum;
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
            ssq += wght[k] * SQUARE(dhat[k]);
        }
        *snew = 0.0;
        for (int k = 0; k < Ndat; k++) {
            dhat[k] *= sqrt(*wsum / ssq);
            *snew += wght[k] * SQUARE(dhat[k] - edis[k]);
        }
        *snew /= *wsum;
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
    }
}
