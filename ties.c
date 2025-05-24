#include "monotone.h"

void primaryApproach(int *ndat, int *blks, double *x, double *w, int *l) {
    int Ndat = *ndat;
    for (int k = 1; k <= Ndat; k++) {
        int blksize = blks[VINDEX(k)];
        if (blksize > 0) {
          double *extracx = (double *)calloc(blksize, sizeof(double)); 
          double *extracw = (double *)calloc(blksize, sizeof(double)); 
          int *extraci = (int *)calloc(blksize, sizeof(int)); 
          for (int i = 1; i <= blksize; i++) {
            extracx[VINDEX(i)] = x[VINDEX(k + i - 1)];
            extracw[VINDEX(i)] = w[VINDEX(k + i - 1)];
            extraci[VINDEX(i)] = l[VINDEX(k + i - 1)];
          }
          (void) mySort(extracx, extraci, &blksize);
          for (int i = 1; i <= blksize; i++) {
            x[VINDEX(k + i - 1)] = extracx[VINDEX(i)];
            w[VINDEX(k + i - 1)] = extracw[VINDEX(i)];
            l[VINDEX(k + i - 1)] = k + extraci[VINDEX(i)] - 1;
          }
          /*
          for (int j = 1; j <= blksize; j++) {
            printf("%d %f %f\n", extraci[VINDEX(j)], extracx[VINDEX(j)], extracw[VINDEX(j)]);
          }
          printf("\n\n");
          */
          free(extracx);
          free(extracw);
          free(extraci);
        }
    }
    (void)monotone(ndat, x, w);
    return;
}

void secondaryApproach(int *ndat, int *blks, double *x, double *w) {
    int nblk = 0, Ndat = *ndat;
    for (int k = 1; k <= Ndat; k++) {
        if (blks[VINDEX(k)] > 0) {
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
    for (int k = 1; k <= nblk; k++) {
        for (int i = l; i <= l + csum[VINDEX(k)] - 1; i++) {
            x[VINDEX(i)] = xave[VINDEX(k)];
        }
        l += csum[VINDEX(k)];
    }
    free(xsum);
    free(wsum);
    free(csum);
    free(xave);
    return;
}

void tertiaryApproach(int *ndat, int *blks, double *x, double *w) {
    int nblk = 0, Ndat = *ndat;
    for (int k = 1; k <= Ndat; k++) {
        if (blks[VINDEX(k)] > 0) {
            nblk++;
        }
    }
    double *xsum = (double *)calloc((size_t)nblk, sizeof(double));
    double *wsum = (double *)calloc((size_t)nblk, sizeof(double));
    double *xave = (double *)calloc((size_t)nblk, sizeof(double));
    double *yave = (double *)calloc((size_t)nblk, sizeof(double));
    int *csum = (int *)calloc((size_t)nblk, sizeof(int));
    for (int k = 1; k <= nblk; k++) {
        yave[VINDEX(k)] = xave[VINDEX(k)];
    }
    (void)tieBlockAverages(ndat, &nblk, blks, x, w, xsum, wsum, csum, xave);
    (void)monotone(&nblk, xave, wsum);
    int l = 1;
    for (int k = 1; k <= nblk; k++) {
        for (int i = l; i <= l + csum[VINDEX(k)] - 1; i++) {
            x[VINDEX(i)] = x[VINDEX(i)] + (xave[VINDEX(k)] - yave[VINDEX(k)]);
        }
        l += csum[VINDEX(k)];
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
    int iblk = 1, Ndat = *ndat, Nblk = *nblk;
    for (int k = 1; k <= Ndat; k++) {
        if (blks[VINDEX(k)] > 0) {
            double sum1 = 0.0, sum2 = 0.0;
            for (int l = k; l <= k + blks[VINDEX(k)] - 1; l++) {
                sum1 += x[VINDEX(l)];
                sum2 += w[VINDEX(l)];
            }
            xsum[VINDEX(iblk)] = sum1;
            wsum[VINDEX(iblk)] = sum2;
            csum[VINDEX(iblk)] = blks[VINDEX(k)];
            iblk++;
        }
    }
    for (int i = 1; i <= Nblk; i++) {
        xave[VINDEX(i)] = xsum[VINDEX(i)] / wsum[VINDEX(i)];
    }
}

int myComp(const void *px, const void *py) {
    double x = ((struct couple *)px)->value;
    double y = ((struct couple *)py)->value;
    return (int)copysign(1.0, x - y);
}

void mySort(double *x, int *k, int *n) {
    int nn = *n;
    struct couple *xi =
        (struct couple *)calloc((size_t)nn, (size_t)sizeof(struct couple));
    for (int i = 0; i < nn; i++) {
        xi[i].value = x[i];
        xi[i].index = i + 1;
    }
    (void)qsort(xi, (size_t)nn, (size_t)sizeof(struct couple), myComp);
    for (int i = 0; i < nn; i++) {
        x[i] = xi[i].value;
        k[i] = xi[i].index;
    }
    free(xi);
}

int main(void) {
	int ndat = 10;
	int blks[10] = {4, 0, 0, 0, 2, 0, 2, 0, 1, 1};
	double x[10] = {4, 3, 2, 1, 3, 1, 3, 6, 4, 3};
	double w[10] = {1, 2, 2, 1, 1, 4, 3, 1, 4, 1};
	int l[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	(void) primaryApproach(&ndat, blks, x, w, l);
	for (int i = 1; i <= ndat; i++) {
	  printf("%f %f %d\n", x[VINDEX(i)], w[VINDEX(i)], l[VINDEX(i)]);
	}
	return EXIT_SUCCESS;
}
