#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define SQUARE(x) ((x) * (x))

void smacofMPInverseV(int* ndat, int* nobj, int* iind, int* jind, double* wght,
                      double* vinv);
void smacofSSMajorize(int* nobj, int* ndim, int* ndat, double* snew, int* iind,
                      int* jind, bool* weighted, double* wght, double* vinv,
                      double* edis, double* dhat, double* xold, double* xnew);

void smacofSSEngine(int* nobj, int* ndim, int* ndat, int* itel, int* ties,
                    int* itmax, int* digits, int* width, bool* verbose,
                    bool* ordinal, bool* weighted, double* sold, double* snew,
                    double* eps, int* iind, int* jind, int* blks, double* wght,
                    double* edis, double* dhat, double* xold, double* xnew);

int main(void) {
    int nobj = 6, ndim = 2, ndat = 15, ties = 1, itmax = 100, digits = 4,
        itel = 1, width = 6;
    bool verbose = true, ordinal = true, weighted = true;
    double sold = 0.0, snew = 0.0, eps = 1e-10;
    int iind[15] = {1, 2, 3, 4, 5, 2, 3, 4, 5, 3, 4, 5, 4, 5, 5};
    int jind[15] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4};
    int blks[15] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double wght[15] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double edis[15] = {2.0, 1.0, 2.0, 2.8, 2.2, 1.0, 2.8, 2.0,
                       2.2, 2.2, 2.2, 2.0, 2.0, 2.0, 1.0};
    double dhat[15] = {1.0, 2.0,  3.0,  4.0,  5.0,  6.0,  7.0, 8.0,
                       9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0};
    double xold[12] = {1.0, -1.0, 0.0, 1.0,  -1.0, 0.0,
                       1.0, 1.0,  1.0, -1.0, -1.0, -1.0};
    double xnew[12] = {0};
    double vinv[15] = {0};
    double sum = 0.0, sam = 0.0, sim = 0.0, wsum = 0.0;
    for (int k = 0; k < ndat; k++) {
        sum += wght[k] * SQUARE(dhat[k]);
        sam += wght[k] * SQUARE(edis[k]);
        wsum += wght[k];
    }
    for (int k = 0; k < ndat; k++) {
        dhat[k] *= sqrt(ndat / sum);
        sim += wght[k] * dhat[k] * edis[k];
    }
    double lbd = sim / sam;
    sum = 0.0;
    for (int k = 0; k < ndat; k++) {
        edis[k] *= lbd;
        sum += wght[k] * SQUARE(dhat[k] - edis[k]);
    }
    for (int k = 0; k < nobj * ndim; k++) {
        xold[k] *= lbd;
    }
    sold = sum / wsum;

    (void)smacofSSEngine(&nobj, &ndim, &ndat, &itel, &ties, 
        &itmax, &digits, &width, &verbose, &ordinal, &weighted,
        &sold, &snew, &eps, iind, jind, blks, wght, edis, dhat, xold, xnew);
    return EXIT_SUCCESS;
}

