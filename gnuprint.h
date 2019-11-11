#ifndef __GNUPRINT_H__
#define __GNUPRINT_H__

#include <stdio.h>
#include "linSystem.h"

int gnuplot(real_t *sol, int nx, int ny, double avg_time, real_t *norm, FILE *fp);

#endif
