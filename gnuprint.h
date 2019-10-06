#ifndef __GNUPRINT_H__
#define __GNUPRINT_H__

#include <stdio.h>
#include "linSystem.h"

int gnuplot(real_t *sol, unsigned int n, double avg_time, real_t *norm, int int_count, FILE *fp);

#endif
