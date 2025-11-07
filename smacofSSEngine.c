#include "smacofSS.h"

void smacofSSEngine(int* nobj, int* ndim, int* ndat, int* itel, int* ties,
                    int* itmax, int* digits, int* width, bool* verbose,
                    bool* ordinal, bool* weighted, double* sold, double* snew,
                    double* eps, int* iind, int* jind, int* blks, double* wght,
                    double* edis, double* dhat, double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    double* vinv = xmalloc(Nobj * (Nobj - 1) * sizeof(double) / 2);
    (void)smacofMPInverseV(nobj, ndat, iind, jind, wght, vinv);
    while (true) {
        (void)smacofSSMajorize(nobj, ndim, ndat, snew, iind, jind, weighted,
                               wght, vinv, edis, dhat, xold, xnew);
        double smid = *snew;
        if (*ordinal) {
            dhat = memcpy(dhat, edis, Ndat * sizeof(double));
            (void)smacofSSMonotone(ndat, ties, snew, iind, jind, blks, edis,
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
}
