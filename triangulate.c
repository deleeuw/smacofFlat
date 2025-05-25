#include <stdlib.h>
#include <stdio.h>

int VINDEX(const int i);
int MINDEX(const int i, const int j, const int n);
int SINDEX(const int i, const int j, const int n);
int TINDEX(const int i, const int j, const int n);

inline int VINDEX(const int i) { return i - 1; }

inline int MINDEX(const int i, const int j, const int n) {
  return (i - 1) + (j - 1) * n;
}

inline int SINDEX(const int i, const int j, const int n) {
  return ((j - 1) * n) - (j * (j - 1) / 2) + (i - j) - 1;
}

inline int TINDEX(const int i, const int j, const int n) {
  return ((j - 1) * n) - ((j - 1) * (j - 2) / 2) + (i - (j - 1)) - 1;
}

void triangulate(double *delta, int *nobj) {
  return;
}

/*
int main(void) {
  int n = 10;
  for (int j = 1; j <= (n - 1); j++) {
    for (int i = (j + 1); i <= n; i++) {
      int k = SINDEX(i, j, n);
      printf("%d %d %d\n", i, j, k);
    }
  }
  return EXIT_SUCCESS;
}
*/