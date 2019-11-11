/*! \file parse.c
    \brief Implementation of argument parsing functions
*/
#include <string.h>
#include <getopt.h>
#include <likwid.h>
#include "linSystem.h"
#include "parse.h"

/*! Structure holds long (>1 char) argument options */
const struct option longopts[] = {
    {"nx", required_argument, 0, 'x'},
    {"ny", required_argument, 0, 'y'}
};

/*!
  \brief Function parses arguments given to main

  \param argc from Main
  \param argv from Main
  \param[out] o structure to that holds the info passed as arguments to main
  \return 0 if parsing was sucessful, 1 otherwise
*/
int parseMain(int argc, char *const argv[], OPT_ARGS_t *o) {
	opterr = 0;
	int c, index = 0;

	while((c = getopt_long_only(argc, argv, "i:o:x:y:", longopts, &index)) != -1) {
		switch(c) {
			case 'i':
				o->MAXITER = atoi(optarg);
				break;
			case 'o':
				strcpy(o->OUTPUT_FILE , optarg);
				break;
            case 'x':
                o->NX = atoi(optarg);
                if(o->NX < 0)
                    return 2;
                break;
            case 'y':
                o->NY = atoi(optarg);
                if(o->NY < 0)
                    return 2;
                break;
			case '?':
				if (optopt == 'i' || optopt == 'o' || strcmp(optarg, "nx") == 0 || strcmp(optarg, "ny") == 0)
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				else
					fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				return 1;
			default:
				return 2;
		}
	}
	return 0;
}

/*!
  \brief Initializes the options struture with default values

  \param opt structure to that holds the indo passed as arguments to main
*/
void initOPTS(OPT_ARGS_t *opt) {
	strcpy(opt->OUTPUT_FILE, "");
    opt->MAXITER = MAXIT;
    opt->NX = 8;
    opt->NY = 8;
}
