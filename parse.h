#include <getopt.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#define MAX_FILE_NAME 15

// Estrutura para tratar a entrada
typedef struct {
	char 	  OUTPUT_FILE[MAX_FILE_NAME];
	int       NX;
    int       NY;
    int       MAXITER;
} OPT_ARGS_t;

int parseMain(int argc, char *const argv[], OPT_ARGS_t *o);
void initOPTS(OPT_ARGS_t *opt);
