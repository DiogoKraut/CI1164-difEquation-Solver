#ifndef __linSystem_H__
#define __linSystem_H__

#define MAXIT 100
#define EPS 1.0e-4

typedef struct {
    float *A;
    float *b;
    unsigned int n;
} linSystem_t;

int gaussSeidel(linSystem_t *S, float *x, float *error);

#endif
