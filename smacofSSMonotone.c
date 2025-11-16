#include "smacofSS.h"

void smacofSSMonotone(int* ndat, int* ties, double *snew, int* iind, int* jind, int* blks,
                      double* edis, double* dhat, double* wght) {
    int Ndat = *ndat;
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
        dhat[k] *= sqrt(((double)Ndat) / ssq);
        *snew += wght[k] * SQUARE(dhat[k] - edis[k]);
    }
    *snew /= Ndat;
}