#include "utils.h"

// Retorna tempo em milisegundos
double timestamp(void)
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}
// 
// void printSys(linSystem_t *S, int nx, int ny) {
//     int i,j;
//     for(i = 0; i < S->n; i++) {
//         for(j = 0; j < S->n; j++) {
//             fprintf(stdout, "% .8f ", S->A[i][j]);
//         }
//         fprintf(stdout, "= % .8f\n", S->b[i]);
//     }
// }
