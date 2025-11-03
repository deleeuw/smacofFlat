#include "smacofSS.h"

void primaryApproach(const size_t *ndat, const size_t *blks, double *x,
                     double *w, double *d, size_t *iind, size_t *jind) {
  size_t Ndat = *ndat;
  for (size_t i = 0; i < Ndat; i++) {
    size_t blksize = blks[i];
    if (blksize > 0) {
      double *extracx = xmalloc(blksize * sizeof(double));
      double *extracw = xmalloc(blksize * sizeof(double));
      double *extracd = xmalloc(blksize * sizeof(double));
      size_t *extraci = xmalloc(blksize * sizeof(size_t));
      size_t *extracj = xmalloc(blksize * sizeof(size_t));
      for (size_t j = 0; j < blksize; j++) {
        extracx[j] = x[i + j];
        extracw[j] = w[i + j];
        extracd[j] = d[i + j];
        extraci[j] = iind[i + j];
        extracj[j] = jind[i + j];
      }
      (void)mySort(extracx, extracw, extracd, extraci, extracj, &blksize);
      for (size_t j = 0; j < blksize; j++) {
        x[i + j] = extracx[j];
        w[i + j] = extracw[j];
        d[i + j] = extracd[j];
        iind[i + j] = extraci[j];
        jind[i + j] = extracj[j];
      }
      xfree(extracx);
      xfree(extracw);
      xfree(extracd);
      xfree(extraci);
      xfree(extracj);
    }
  }
  double *ww = xmalloc(Ndat * sizeof(double));
  for (size_t i = 0; i < Ndat; i++) {
    ww[i] = w[i];
  }
  (void)monotone(ndat, x, ww);
  xfree(ww);
  return;
}

void secondaryApproach(const size_t *ndat, const size_t *blks, double *x,
                       double *w) {
  size_t nblk = 0, Ndat = *ndat;
  for (size_t k = 0; k < Ndat; k++) {
    if (blks[k] > 0) {
      nblk++;
    }
  }
  double *xsum = xmalloc(nblk * sizeof(double));
  double *wsum = xmalloc(nblk * sizeof(double));
  double *xave = xmalloc(nblk * sizeof(double));
  size_t *csum = xmalloc(nblk * sizeof(size_t));
  (void)tieBlockAverages(ndat, &nblk, blks, x, w, xsum, wsum, csum, xave);
  (void)monotone(&nblk, xave, wsum);
  size_t l = 1;
  for (size_t k = 0; k < nblk; k++) {
    for (size_t i = l; i <= l + csum[k] - 1; i++) {
      x[i - 1] = xave[k];
    }
    l += csum[k];
  }
  xfree(xsum);
  xfree(wsum);
  xfree(csum);
  xfree(xave);
  return;
}

void tertiaryApproach(const size_t *ndat, const size_t *blks, double *x,
                      double *w) {
  size_t nblk = 0, Ndat = *ndat;
  for (size_t k = 0; k < Ndat; k++) {
    if (blks[k] > 0) {
      nblk++;
    }
  }
  double *xsum = xmalloc(nblk * sizeof(double));
  double *wsum = xmalloc(nblk * sizeof(double));
  double *xave = xmalloc(nblk * sizeof(double));
  double *yave = xmalloc(nblk * sizeof(double));
  size_t *csum = xmalloc(nblk * sizeof(size_t));
  (void)tieBlockAverages(ndat, &nblk, blks, x, w, xsum, wsum, csum, xave);
  for (size_t k = 0; k < nblk; k++) {
    yave[k] = xave[k];
  }
  (void)monotone(&nblk, xave, wsum);
  size_t l = 1;
  for (size_t k = 0; k < nblk; k++) {
    for (size_t i = l; i <= l + csum[k] - 1; i++) {
      x[i - 1] = xave[k] + (x[i - 1] - yave[k]);
    }
    l += csum[k];
  }
  xfree(xsum);
  xfree(wsum);
  xfree(xave);
  xfree(yave);
  xfree(csum);
  return;
}

void tieBlockAverages(const size_t *ndat, const size_t *nblk,
                      const size_t *blks, const double *x, const double *w,
                      double *xsum, double *wsum, size_t *csum, double *xave) {
  size_t iblk = 0, Ndat = *ndat, Nblk = *nblk;
  for (size_t k = 0; k < Ndat; k++) {
    if (blks[k] > 0) {
      double sum1 = 0.0, sum2 = 0.0;
      for (size_t l = k; l < k + blks[k]; l++) {
        sum1 += w[l] * x[l];
        sum2 += w[l];
      }
      xsum[iblk] = sum1;
      wsum[iblk] = sum2;
      csum[iblk] = blks[k];
      iblk++;
    }
  }
  for (size_t i = 0; i < Nblk; i++) {
    xave[i] = xsum[i] / wsum[i];
  }
}

void monotone(const size_t *n, double *x, double *w) {
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
    for (size_t i = from; i >= to; i--)
      rx[i] = xk;
    from = to - 1;
  }
  xfree(idx);
}
