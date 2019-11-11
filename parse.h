#include <getopt.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#define MAX_FILE 15

/*! Struture for handling arguments */
typedef struct {
	char 	  OUTPUT_FILE[MAX_FILE]; /*!< Specified output file. Stdout if empty    */
	int       NX;                    /*!< Number of subdivisions in the X dimension */
  int       NY;                    /*!< Number of subdivisions in the Y dimension */
  int       MAXITER;               /*!< Guass-Seidel iteration limit              */
} OPT_ARGS_t;

int parseMain(int argc, char *const argv[], OPT_ARGS_t *o);
void initOPTS(OPT_ARGS_t *opt);
