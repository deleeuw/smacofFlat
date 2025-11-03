#ifndef SMACOF_SS_H
#define SMACOF_SS_H

#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SQUARE(x) ((x) * (x))

void primaryApproach(const size_t *, const size_t *, double *, double *,
                     double *, size_t *, size_t *);
void secondaryApproach(const size_t *, const size_t *, double *, double *);
void tertiaryApproach(const size_t *, const size_t *, double *, double *);
void tieBlockAverages(const size_t *, const size_t *, const size_t *,
                      const double *, const double *, double *, double *,
                      size_t *, double *);
void monotone(const size_t *, double *, double *);
int myComp(const void *, const void *);
void mySort(double *, double *, double *, size_t *, size_t *, const size_t *);

static inline void *xmalloc(const size_t size) {
  void *p = malloc(size);
  if (!p && size) {
    fprintf(stderr, "FATAL: malloc(%zu) failed\n", size);
    abort();
  }
  return p;
}

static inline void *xcalloc(const size_t nmemb, const size_t size) {
  /* basic overflow guard */
  if (size && nmemb > SIZE_MAX / size) {
    fprintf(stderr, "FATAL: calloc overflow (%zu,%zu)\n", nmemb, size);
    abort();
  }
  void *p = calloc(nmemb, size);
  if (!p && nmemb && size) {
    fprintf(stderr, "FATAL: calloc(%zu,%zu) failed\n", nmemb, size);
    abort();
  }
  return p;
}

static inline void *xrealloc(void *ptr, const size_t size) {
  void *p = realloc(ptr, size);
  if (!p && size != 0) {
    fprintf(stderr, "FATAL: realloc(%p,%zu) failed\n", ptr, size);
    abort();
  }
  return p;
}

#define xfree(p)                                                               \
  {                                                                            \
    if ((p) != NULL) {                                                         \
      free(p);                                                                 \
      p = NULL;                                                                \
    }                                                                          \
  }

#endif /* SMACOF_SS_H */