#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define VINDEX(i) (i - 1)

void primaryApproach(int *, int *, double *, double *, int *); 
void secondaryApproach(int *, int *, double *, double *);
void tertiaryApproach(int *, int *, double *, double *);
void tieBlockAverages(int *, int *, int *, double *, double *,
                      double *, double *, int *, double *);
void monotone(int*, double*, double*);
int myComp(const void *, const void *);
void mySort(double *, int *, int *);

struct couple {
    double value;
    int index;
};