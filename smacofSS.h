#ifndef SMACOF_SS_H
#define SMACOF_SS_H

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>

#define SQUARE(x) ((x) * (x))

void primaryApproach(int *, int *, double *, double *, double *, int *, int *);
void secondaryApproach(int *, int *, double *, double *);
void tertiaryApproach(int *, int *, double *, double *);
void tieBlockAverages(int *, int *, int *, double *, double *, double *,
                      double *, int *, double *);
void monotone(int *, double *, double *);
int myComp(const void *, const void *);
void mySort(double *, double *, double *, int *, int *, int *);

#endif /* SMACOF_SS_H */