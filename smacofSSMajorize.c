#include "smacofSS.h"

void smacofSSMajorize(const int* nobj, const int* ndim, const int* ndat,
                      const int* itel, int* kord, const int* nord, int* iind,
                      int* jind, const int* iord, const int* safe,
                      const int* weighted, double* wght, double* vinv,
                      double* dhat, double* xold, double* xnew) {
    int Nobj = *nobj, Ndim = *ndim, Itel = *itel, Nord = *nord, Ndat = *ndat,
        Kord = *kord;
    Kord = iord[(Itel - 1) % Nord];
    int ncor = Nobj * Ndim, Safe = *safe;
    double snul = 0.0, snew = 0.0;
    double* xone = xmalloc(Nobj * Ndim * sizeof(double));
    double* xtwo = xmalloc(Nobj * Ndim * sizeof(double));
    double* xthr = xmalloc(Nobj * Ndim * sizeof(double));
    double* done = xmalloc(Nobj * Ndim * sizeof(double));
    double* dtwo = xmalloc(Nobj * Ndim * sizeof(double));
    double* dthr = xmalloc(Nobj * Ndim * sizeof(double));
    double* edis = xmalloc(Ndat * sizeof(double));
    (void)smacofSSGuttmanTransform(nobj, ndim, ndat, iind, jind, weighted, wght,
                                   vinv, dhat, xold, xone);
    if (Safe) {
        (void)smacofSSDistances(nobj, ndim, ndat, iind, jind, xone, edis);
        snul = smacofSSLoss(ndat, edis, dhat, wght);
    }
    if (Kord == 0) {
        for (int i = 0; i < ncor; i++) {
            xnew[i] = xone[i];
        }
    } else {
        for (int i = 0; i < ncor; i++) {
            done[i] = xone[i] - xold[i];
        }
        if (Kord == 1) {
            for (int i = 0; i < ncor; i++) {
                xnew[i] = xold[i] + 2.0 * done[i];
            }
        } else {
            (void)smacofSSGuttmanTransform(nobj, ndim, ndat, iind, jind,
                                           weighted, wght, vinv, dhat, xone,
                                           xtwo);
            for (int i = 0; i < ncor; i++) {
                dtwo[i] = xtwo[i] - 2.0 * xone[i] + xold[i];
            }
            if (Kord == 2) {
                double sum1 = 0.0, sum2 = 0.0;
                for (int i = 0; i < ncor; i++) {
                    sum1 += done[i] * dtwo[i];
                    sum2 += SQUARE(dtwo[i]);
                }
                double sig2 = fabs(sum1 / sum2);
                for (int i = 0; i < ncor; i++) {
                    xnew[i] =
                        xold[i] + 2 * sig2 * done[i] + SQUARE(sig2) * dtwo[i];
                }
            } else {
                (void)smacofSSGuttmanTransform(nobj, ndim, ndat, iind, jind,
                                               weighted, wght, vinv, dhat, xtwo,
                                               xthr);
                for (int i = 0; i < ncor; i++) {
                    dthr[i] = xthr[i] - 3.0 * xtwo[i] + 3.0 * xone[i] - xold[i];
                }
                if (Kord == 3) {
                    double sum1 = 0.0, sum2 = 0.0;
                    for (int i = 0; i < ncor; i++) {
                        sum1 += dtwo[i] * dthr[i];
                        sum2 += SQUARE(dthr[i]);
                    }
                    double sig3 = fabs(sum1 / sum2);
                    for (int i = 0; i < ncor; i++) {
                        xnew[i] = xold[i] + 3.0 * sig3 * done[i] +
                                  3.0 * SQUARE(sig3) * dtwo[i] +
                                  CUBE(sig3) * dthr[i];
                    }
                }
            }
        }
    }
    if (Safe) {
        (void)smacofSSDistances(nobj, ndim, ndat, iind, jind, xnew, edis);
        snew = smacofSSLoss(ndat, edis, dhat, wght);
        if (snul < snew) {
            for (int i = 0; i < ncor; i++) {
                xnew[i] = xone[i];
                Kord = 0;
            }
        }
    }
    xfree(xone);
    xfree(xtwo);
    xfree(xthr);
    xfree(done);
    xfree(dtwo);
    xfree(dthr);
    xfree(edis);
    *kord = Kord;
    return;
}

void smacofSSGuttmanTransform(const int* nobj, const int* ndim, const int* ndat,
                              int* iind, int* jind, const int* weighted,
                              double* wght, double* vinv, double* dhat,
                              double* xold, double* xnew) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim, Weighted = *weighted;
    double* xtmp = xmalloc(Nobj * Ndim * sizeof(double));
    double* edis = xmalloc(Ndat * sizeof(double));
    (void)smacofSSDistances(nobj, ndim, ndat, iind, jind, xold, edis);
    for (int k = 0; k < Nobj * Ndim; k++) {
        xtmp[k] = 0.0;
    }
    for (int k = 0; k < Ndat; k++) {
        int i = iind[k], j = jind[k];
        if (edis[k] < EPS) {
          continue;
        }
        double ecof = wght[k] * dhat[k] / edis[k];
        for (int s = 0; s < Ndim; s++) {
            int iobj = i + Nobj * s, jobj = j + Nobj * s;
            double add = ecof * (xold[iobj] - xold[jobj]);
            xtmp[iobj] += add;
            xtmp[jobj] -= add;
        }
    }
    for (int k = 0; k < Nobj * Ndim; k++) {
        if (Weighted) {
            xnew[k] = 0.0;
        } else {
            xnew[k] = xtmp[k];
        }
    }
    // here follows a bug
    if (Weighted) {
        /*
         for (int k = 0; k < Ndat; k++) {
            int i = iind[k], j = jind[k];
            double ecof = -vinv[k];
            for (int s = 0; s < Ndim; s++) {
                int iobj = i + Nobj * s, jobj = j + Nobj * s;
                double add = ecof * (xtmp[iobj] - xtmp[jobj]);
                xnew[iobj] += add;
                xnew[jobj] -= add;
            }
        }
        */
      int k = 0;
      for (int j = 0; j < Nobj - 1; j++) {
        for (int i = j + 1; i < Nobj; i++) {
          double ecof = -vinv[k];
          for (int s = 0; s < Ndim; s++) {
            int iobj = i + Nobj * s, jobj = j + Nobj * s;
            double add = ecof * (xtmp[iobj] - xtmp[jobj]);
            xnew[iobj] += add;
            xnew[jobj] -= add;
          }          
          k++;
        }
      }
    } else {
        for (int i = 0; i < Nobj; i++) {
            int is = i;
            for (int s = 0; s < Ndim; s++) {
                xnew[is] /= (double)Nobj;
                is += Nobj;
            }
        }
    }
    xfree(xtmp);
    xfree(edis);
    return;
}

void smacofSSDistances(const int* nobj, const int* ndim, const int* ndat,
                       int* iind, int* jind, double* xmat, double* edis) {
    int Ndat = *ndat, Nobj = *nobj, Ndim = *ndim;
    for (int k = 0; k < Ndat; k++) {
        int is = iind[k], js = jind[k];
        double sum = 0.0;
        for (int s = 0; s < Ndim; s++) {
            sum += SQUARE(xmat[is] - xmat[js]);
            is += Nobj;
            js += Nobj;
        }
        edis[k] = sqrt(sum);
    }
}