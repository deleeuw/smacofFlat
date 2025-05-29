
void primaryApproach(int *ndat, int *blks, double *x, double *w, int *indi) {
    int Ndat = *ndat;
    for (int i = 0; i < Ndat; i++) {
        int blksize = blks[i];
        if (blksize > 0) {
            double *extracx = (double *)calloc(blksize, sizeof(double));
            double *extracw = (double *)calloc(blksize, sizeof(double));
            int *extraci = (int *)calloc(blksize, sizeof(int));
            for (int j = 0; j < blksize; j++) {
                extracx[j] = x[i + j];
                extracw[j] = w[i + j];
                extraci[j] = indi[i + j];
            }
            (void)mySort(extracx, extracw, extraci, &blksize);
            for (int j = 0; j < blksize; j++) {
                x[i + j] = extracx[j];
                w[i + j] = extracw[j];
                indi[i + j] = extraci[j];
            }
            free(extracx);
            free(extracw);
            free(extraci);
        }
    }
    double *ww = (double *)calloc(Ndat, sizeof(double));
    for (int i = 0; i < Ndat; i++) {
        ww[i] = w[i];
    }
    (void)monotone(ndat, x, ww);
    free(ww);
    return;
}

void secondaryApproach(int *ndat, int *blks, double *x, double *w) {
    int nblk = 0, Ndat = *ndat;
    for (int k = 0; k < Ndat; k++) {
        if (blks[k] > 0) {
            nblk++;
        }
    }
    double *xsum = (double *)calloc((size_t)nblk, sizeof(double));
    double *wsum = (double *)calloc((size_t)nblk, sizeof(double));
    double *xave = (double *)calloc((size_t)nblk, sizeof(double));
    int *csum = (int *)calloc((size_t)nblk, sizeof(int));
    (void)tieBlockAverages(ndat, &nblk, blks, x, w, xsum, wsum, csum, xave);
    (void)monotone(&nblk, xave, wsum);
    int l = 1;
    for (int k = 0; k < nblk; k++) {
        for (int i = l; i <= l + csum[k] - 1; i++) {
            x[i - 1] = xave[k];
        }
        l += csum[k];
    }
    free(xsum);
    free(wsum);
    free(csum);
    free(xave);
    return;
}

void tertiaryApproach(int *ndat, int *blks, double *x, double *w) {
    int nblk = 0, Ndat = *ndat;
    for (int k = 0; k < Ndat; k++) {
        if (blks[k] > 0) {
            nblk++;
        }
    }
    double *xsum = (double *)calloc((size_t)nblk, sizeof(double));
    double *wsum = (double *)calloc((size_t)nblk, sizeof(double));
    double *xave = (double *)calloc((size_t)nblk, sizeof(double));
    double *yave = (double *)calloc((size_t)nblk, sizeof(double));
    int *csum = (int *)calloc((size_t)nblk, sizeof(int));
    (void)tieBlockAverages(ndat, &nblk, blks, x, w, xsum, wsum, csum, xave);
    for (int k = 0; k < nblk; k++) {
        yave[k] = xave[k];
    }
    (void)monotone(&nblk, xave, wsum);
    int l = 1;
    for (int k = 0; k < nblk; k++) {
        for (int i = l; i <= l + csum[k] - 1; i++) {
            x[i - 1] = xave[k] + (x[i - 1] - yave[k]);
        }
        l += csum[k];
    }
    free(xsum);
    free(wsum);
    free(xave);
    free(yave);
    free(csum);
    return;
}

void tieBlockAverages(int *ndat, int *nblk, int *blks, double *x, double *w,
                      double *xsum, double *wsum, int *csum, double *xave) {
    int iblk = 0, Ndat = *ndat, Nblk = *nblk;
    for (int k = 0; k < Ndat; k++) {
        if (blks[k] > 0) {
            double sum1 = 0.0, sum2 = 0.0;
            for (int l = k; l < k + blks[k]; l++) {
                sum1 += w[l] * x[l];
                sum2 += w[l];
            }
            xsum[iblk] = sum1;
            wsum[iblk] = sum2;
            csum[iblk] = blks[k];
            iblk++;
        }
    }
    for (int i = 0; i < Nblk; i++) {
        xave[i] = xsum[i] / wsum[i];
    }
}

int myComp(const void *px, const void *py) {
    double x = ((struct triple *)px)->value;
    double y = ((struct triple *)py)->value;
    return (int)copysign(1.0, x - y);
}

void mySort(double *x, double *w, int *indi, int *n) {
    int nn = *n;
    struct triple *xi =
        (struct triple *)calloc((size_t)nn, (size_t)sizeof(struct triple));
    for (int i = 0; i < nn; i++) {
        xi[i].value = x[i];
        xi[i].weight = w[i];
        xi[i].index = indi[i];
    }
    (void)qsort(xi, (size_t)nn, (size_t)sizeof(struct triple), myComp);
    for (int i = 0; i < nn; i++) {
        x[i] = xi[i].value;
        w[i] = xi[i].weight;
        indi[i] = xi[i].index;
    }
    free(xi);
}

#include "monotone.h"

void monotone(int *n, double *x, double *w)
// Function monotone(),
// performs simple linear ordered monotone regression
// Copyright (C) 2020 Frank M.T.A. Busing (e-mail: busing at fsw dot leidenuniv
// dot nl) This function is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation. This program is distributed in the hope that it
// will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
// Public License for more details. You should have received a copy of the GNU
// General Public License along with this function. If not, see
// <https://www.gnu.org/licenses/>.
{
    double *rx = &x[-1];
    double *rw = &w[-1];
    size_t *idx = (size_t *)calloc(*n + 1, sizeof(size_t));
    idx[0] = 0;
    idx[1] = 1;
    size_t b = 1;
    double xbm1 = rx[b];
    double wbm1 = rw[b];
    for (size_t i = 2; i <= *n; i++) {
        b++;
        double xb = rx[i];
        double wb = rw[i];
        if (xbm1 > xb) {
            b--;
            double sb = wbm1 * xbm1 + wb * xb;
            wb += wbm1;
            xb = sb / wb;
            while (i < *n && xb >= rx[i + 1]) {
                i++;
                sb += rw[i] * rx[i];
                wb += rw[i];
                xb = sb / wb;
            }
            while (b > 1 && rx[b - 1] > xb) {
                b--;
                sb += rw[b] * rx[b];
                wb += rw[b];
                xb = sb / wb;
            }
        }
        rx[b] = xbm1 = xb;
        rw[b] = wbm1 = wb;
        idx[b] = i;
    }
    size_t from = *n;
    for (size_t k = b; k > 0; k--) {
        const size_t to = idx[k - 1] + 1;
        const double xk = rx[k];
        for (size_t i = from; i >= to; i--) rx[i] = xk;
        from = to - 1;
    }
    free(idx);
}  // monotone

/*
int main(void) {
    double x[9] = {8.0, 3.0, 5.0, 2.0, 1.0, 9.0, 6.0, 7.0, 4.0};
    double w[9] = {1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 1.0};
    int blks[9] = {3, 0, 0, 3, 0, 0, 2, 0, 1};
    int iind[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    int ndat = 9, nblk = 4;
    for (int i = 0; i < ndat; i++) {
        printf("* %f %f %d %d\n", x[i], w[i], blks[i], iind[i]);
    }
    printf("\n\n");
    double *xsum = (double *)calloc(nblk, sizeof(double));
    double *wsum = (double *)calloc(nblk, sizeof(double));
    double *xave = (double *)calloc(nblk, sizeof(double));
    int *csum = (int *)calloc(nblk, sizeof(int));
    (void)tieBlockAverages(&ndat, &nblk, blks, x, w, xsum, wsum, csum, xave);
    for (int i = 0; i < nblk; i++) {
        printf("%f %f %f %d\n", xsum[i], wsum[i], xave[i], csum[i]);
    }
    printf("\n\n");
    free(xsum);
    free(wsum);
    free(xave);
    free(csum);
    (void)secondaryApproach(&ndat, blks, x, w);
    for (int i = 0; i < ndat; i++) {
        printf("%f %f %d\n", x[i], w[i], iind[i]);
    }
    return EXIT_SUCCESS;
}
*/