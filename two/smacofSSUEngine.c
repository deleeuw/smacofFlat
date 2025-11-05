#include "smacofSS.h"

void smacofSSUEngine(int* nobj, int* ndim, int* ndat, int* itel, int* ties,
                     int* itmax, int* digits, int* width, bool* verbose,
                     bool* ordinal, double* sold, double* snew, double* eps,
                     int* iind, int* jind, int* blks, double* edis,
                     double* dhat, double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    double* wght = xcalloc(Ndat, sizeof(double));
    for (int k = 0; k < Ndat; k++) {
        wght[k] = 1.0;
    }
    while (true) {
        (void)smacofSSUMajorize(nobj, ndim, ndat, snew, iind, jind, edis, dhat,
                                xold, xnew);
        double smid = *snew;
        if (*ordinal) {
            dhat = memcpy(dhat, edis, Ndat * sizeof(double));
            (void)smacofSSUMonotone(ndat, ties, snew, iind, jind, blks, edis,
                                    dhat, wght);
        }
        if (*verbose) {
            if (*ordinal) {
                printf("itel %4d sold %*.*f smid %*.*f snew %*.*f\n", *itel,
                       *width, *digits, *sold, *width, *digits, smid, *width,
                       *digits, *snew);
            } else {
                printf("itel %4d sold %*.*f snew %*.*f\n", *itel, *width,
                       *digits, *sold, *width, *digits, *snew);
            }
        }
        if ((*itel == *itmax) || ((*sold - *snew) < *eps)) {
            break;
        }
        xold = memcpy(xold, xnew, Nobj * Ndim * sizeof(double));
        *sold = *snew;
        *itel += 1;
    }
    xfree(wght);
}
