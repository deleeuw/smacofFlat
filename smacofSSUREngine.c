#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#define SQUARE(x) ((x) * (x))

void smacofStepUnweighted(int *, int *, int *, int *, int *, double *, double *, double *, double *);
                                                        
void smacofSSUREngine(int *nobj,
                      int *ndim,
                      int *ndat,
                      int *itel, 
                      int *itmax,
                      bool *verbose,
                      double *sold, 
                      double *snew,
                      double *eps,
                      int *iind,
                      int *jind,
                      double *edis,
                      double *dhat,
                      double *xold,
                      double *xnew){
  for (;;) {
    for (int i = 0; i < *ndat; i++) {
      printf(" %f ", edis[i]);
    }
    printf("\n");
    (void) smacofStepUnweighted(nobj, ndim, ndat, iind, jind, edis, dhat, xold, xnew);
    *snew = 0.0;
    for (int i = 0; i < *ndat; i++) {
      *snew += SQUARE(dhat[i] - edis[i]);
      printf(" %f ", edis[i]);
    }
    printf("\n");
    if (*verbose) {
      printf("itel %4d sold %15.10f snew %15.10f\n", *itel, *sold, *snew);
      }
    if ((*itel == *itmax) || ((*sold - *snew) < *eps)) {
      break;
    }
    *sold = *snew;
    *itel += 1;
  }
}