#ifndef __linSystem_H__
#define __linSystem_H__

#define MAXIT 10000
#define EPS 1.0e-4

typedef float real_t; // typedef for switching between float/double

typedef struct {
    real_t **A;
    real_t *b;
    unsigned int n;
} linSystem_t;

int gaussSeidel(linSystem_t *S, real_t *x, real_t error);

#endif
