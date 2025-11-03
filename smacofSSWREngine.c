#include "smacofSS.h"

void smacofSSWREngine(const size_t* nobj, const size_t* ndim, const size_t* ndat, size_t* itel, const size_t* itmax,
                      const int* digits, const int* width, const bool* verbose, double* wsum,
                      double* sold, double* snew, const double* eps, const size_t* iind,
                      const size_t* jind, double* edis, const double* dhat, const double* wght,
                      const double* vinv, double* xold, double* xnew) {
  size_t Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
  while (true) {
    double* xtmp = xmalloc(Nobj * Ndim * sizeof(double));
    for (size_t k = 0; k < Ndat; k++) {
      if (edis[k] == 0.0) {
        continue;
      }
      size_t is = iind[k] - 1, js = jind[k] - 1;
      double elem = wght[k] * dhat[k] / edis[k];
      for (size_t s = 0; s < Ndim; s++) {
        double add = elem * (xold[is] - xold[js]);
        xtmp[is] += add;
        xtmp[js] -= add;
        is += Nobj;
        js += Nobj;
      }
    }
    for (size_t k = 0; k < Nobj * Ndim; k++) {
      xnew[k] = 0.0;
    }
    size_t k = 0;
    for (size_t j = 0; j < Nobj - 1; j++) {
      for (size_t i = j + 1; i < Nobj; i++) {
        double elem = vinv[k];
        size_t is = i, js = j;
        for (size_t s = 0; s < Ndim; s++) {
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
    for (size_t k = 0; k < Ndat; k++) {
      size_t is = iind[k] - 1, js = jind[k] - 1;
      double sum = 0.0;
      for (size_t s = 0; s < Ndim; s++) {
        sum += SQUARE(xnew[is] - xnew[js]);
        is += Nobj;
        js += Nobj;
      }
      edis[k] = sqrt(sum);
    }
    *snew = 0.0;
    for (size_t k = 0; k < Ndat; k++) {
      *snew += wght[k] * SQUARE(dhat[k] - edis[k]);
    }
    *snew /= *wsum;
    if (*verbose) {
      printf("itel %4zu sold %*.*f snew %*.*f\n", *itel, *width, *digits, *sold,
             *width, *digits, *snew);
    }
    if ((*itel == *itmax) || ((*sold - *snew) < *eps)) {
      break;
    }
    xold = memcpy(xold, xnew, Nobj * Ndim * sizeof(double));
    *sold = *snew;
    *itel += 1;
  }
}