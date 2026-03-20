#include <stdio.h>
#include <stdlib.h>

size_t  nrow = 10, ncol = 5;

int main(void) {
    double (*a)[ncol] = calloc(nrow * ncol, sizeof(double));
    for (size_t i = 0; i < nrow; i++) {
        for (size_t j = 0; j < ncol; j++) {
            a[i][j] = (i + 1) * (j + 1);
            printf(" %6.0f ", a[i][j]);
        }
        printf("\n");
    }
    free(a);
}