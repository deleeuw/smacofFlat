
#ifndef SMACOF_ALLOC_H
#define SMACOF_ALLOC_H

/* Abort-on-OOM (Out of Memory) policy (common for daemons/servers).
   If you need non-fatal OOM, use *_try variants separately. */

static inline void *xmalloc(size_t size) {
    void *p = malloc(size);
    if (!p && size) {
        fprintf(stderr, "FATAL: malloc(%zu) failed\n", size);
        abort();
    }
    return p;
}

static inline void *xmalloc_try(size_t size) {
    void *p = malloc(size);
    if (!p && size) {
        return NULL;
    }
    return p;
}

static inline void *xcalloc(size_t nmemb, size_t size) {
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

static inline void *xcalloc_try(size_t nmemb, size_t size) {
    /* basic overflow guard */
    if (size && nmemb > SIZE_MAX / size) {
        return NULL;
    }
    void *p = calloc(nmemb, size);
    if (!p && nmemb && size) {
        return NULL;
    }
    return p;
}


static inline void *xrealloc(void *ptr, size_t size) {
    void *p = realloc(ptr, size);
    if (!p && size != 0) {
        fprintf(stderr, "FATAL: realloc(%p,%zu) failed\n", ptr, size);
        abort();
    }
    return p;
}

static inline void *xrealloc_try(void *ptr, size_t size) {
    void *p = realloc(ptr, size);
    if (!p && size != 0) {
        return NULL;
    }
    return p;
}

/* Free + NULL 
*  We put scope in here to avoid dangling if-else on non braces statement
*/
#define xfree(p) { if ((p) != NULL) { free(p); p = NULL; } }

/* Bounded str duplicate: must provide valid cstr with NUL-terminate. */
static inline char *xstrdup(const char *s) {
    if (!s) {
        /* Define our policy: duplicate NULL → empty string */
        char *z = xmalloc(1);
        z[0] = '\0';
        return z;
    }
    size_t n = strlen(s);
    /* +1 checked overflow */
    if (n >= SIZE_MAX) {
        fprintf(stderr, "FATAL: xstrdup overflow\n");
        abort();
    }

    char *p = xmalloc(n + 1);
    memcpy(p, s, n + 1); /* includes '\0' */
    return p;
}

/* Bounded str duplicate: copy at most n bytes, n is size-1 (withouth NUL-terminate). */
static inline char *xstrndup(const char *s, size_t n) {
    if (!s) {
        char *z = xmalloc(1);
        z[0] = '\0';
        return z;
    }

    size_t m = strlen(s);

    // size_t m = strnlen(s, maxlen);
    if (m >= SIZE_MAX) {
        fprintf(stderr, "FATAL: xstrndup overflow\n");
        abort();
    }

    char *p = xmalloc(m + 1);

    if (m){
        memcpy(p, s, m);
    }

    p[m] = '\0';
    return p;
}

/* Binary memcopy helper. */
static inline void *xmemcpy(const void *src, size_t n) {
    if (!src && n) {
        return NULL;
    }

    void *p = xmalloc(n ? n : 1);
    if (n){ 
        memcpy(p, src, n);
    }

    return p;
}

#endif /* SMACOF_ALLOC.H */