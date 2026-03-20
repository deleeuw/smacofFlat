#include "smacofSS.h"

void smacofSSSMatrixPrint(double* mat, int* nobj, int* ndat, int* iind,
                          int* jind, int* width, int* digits) {
    int Nobj = *nobj, Ndat = *ndat;
    double* out = xcalloc(Nobj * Nobj, sizeof(double));
    for (int k = 0; k < Ndat; k++) {
        int i = iind[k];
        int j = jind[k];
        out[i + Nobj * j] = out[j + Nobj * i] = mat[k];
    }
    for (int i = 0; i < Nobj; i++) {
        for (int j = 0; j < Nobj; j++) {
            printf(" %*.*f ", *width, *digits, out[i + Nobj * j]);
        }
        printf("\n");
    }
    printf("\n\n");
    return;
}

void smacofSSRMatrixPrint(double* mat, int* nrow, int* ncol, int* width,
                          int* digits) {
    int Nrow = *nrow, Ncol = *ncol;
    for (int i = 0; i < Nrow; i++) {
        for (int j = 0; j < Ncol; j++) {
            int k = i + j * Nrow;
            printf(" %*.*f ", *width, *digits, mat[k]);
        }
        printf("\n");
    }
    printf("\n\n");
}

void smacofSSVectorPrint(double* vec, int* n, int* width, int* digits) {
    for (int i = 0; i < *n; i++) {
        printf(" %*.*f ", *width, *digits, vec[i]);
    }
    printf("\n\n");
    return;
}