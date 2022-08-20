#include <unistd.h>

#include "Options.h"

void PrintOpts(const OPTS &opts) {

	printf("\nUsage:   ./get_features [-option] [argument]\n\n");
	printf("Options:  -i input.pdb  \\       # (required) input PDB file\n");
	printf("          -r reference.pdb \\    # (optional) reference PDB file\n");
	printf("          -j output.json \\      # (optional) atomic features\n");
	printf("          -p output.pdb \\       # (optional) cleaned model\n");
	printf("          -d DMAX \\             # (%.2lf) contact distance\n", opts.dmax);
	printf("          -t TOPN \\             # (%d) top contacts to save\n", opts.topn);
	printf("          -v VERBOSITY \\        # (%d) verbosity level\n", opts.verb);
	printf("\n");

}

bool GetOpts(int argc, char *argv[], OPTS &opts) {

	char tmp;
	while ((tmp = getopt(argc, argv, "hi:r:j:d:p:v:t:")) != -1) {
		switch (tmp) {
		case 'h': /* help */
			printf("!!! HELP !!!\n");
			return false;
			break;
		case 'i': /* model PDB*/
			opts.model = std::string(optarg);
			break;
		case 'r': /* reference PDB */
			opts.ref = std::string(optarg);
			break;
		case 'j': /* output file */
			opts.json = std::string(optarg);
			break;
		case 'd': /* chain ID */
			opts.dmax = atof(optarg);
			break;
		case 't':
			opts.topn = atoi(optarg);
			break;
		case 'p':
			opts.pdb = std::string(optarg);
			break;
		case 'v':
			opts.verb = atoi(optarg);
			break;
		default:
			return false;
			break;
		}
	}

	if (opts.model == "") {
		printf("Error: input file not specified ('-i')\n");
		return false;
	}

	return true;

}

