#include "smacofSS.h"
#include "smacofAlloc.h"

void primaryApproach(int *ndat, int *blks, double *x, double *w, double *d,
                     int *iind, int *jind) {
    int Ndat = *ndat;
    for (int i = 0; i < Ndat; i++) {
        int blksize = blks[i];
        if (blksize > 0) {
            double *extracx = xcalloc(blksize, sizeof(double));
            double *extracw = xcalloc(blksize, sizeof(double));
            double *extracd = xcalloc(blksize, sizeof(double));
            int *extraci = xcalloc(blksize, sizeof(int));
            int *extracj = xcalloc(blksize, sizeof(int));
            for (int j = 0; j < blksize; j++) {
                extracx[j] = x[i + j];
                extracw[j] = w[i + j];
                extracd[j] = d[i + j];
                extraci[j] = iind[i + j];
                extracj[j] = jind[i + j];
            }
            (void)mySort(extracx, extracw, extracd, extraci, extracj, &blksize);
            for (int j = 0; j < blksize; j++) {
                x[i + j] = extracx[j];
                w[i + j] = extracw[j];
                d[i + j] = extracd[j];
                iind[i + j] = extraci[j];
                jind[i + j] = extracj[j];
            }
            free(extracx);
            free(extracw);
            free(extracd);
            free(extraci);
            free(extracj);
        }
    }
    double *ww = xcalloc(Ndat, sizeof(double));
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
    double *xsum = xcalloc((size_t)nblk, sizeof(double));
    double *wsum = xcalloc((size_t)nblk, sizeof(double));
    double *xave = xcalloc((size_t)nblk, sizeof(double));
    int *csum = xcalloc((size_t)nblk, sizeof(int));
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
    double *xsum = xcalloc((size_t)nblk, sizeof(double));
    double *wsum = xcalloc((size_t)nblk, sizeof(double));
    double *xave = xcalloc((size_t)nblk, sizeof(double));
    double *yave = xcalloc((size_t)nblk, sizeof(double));
    int *csum = xcalloc((size_t)nblk, sizeof(int));
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

void monotone(int *n, double *x, double *w)
{
    double *rx = &x[-1];
    double *rw = &w[-1];
    size_t *idx = xcalloc(*n + 1, sizeof(size_t));
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
}  
