#include <stdlib.h>

#include "export.h"
#include "forward.h"
#include "init.h"
#include "memory.h"
#include "parse_args.hpp"
#include "shalw.h"

double *	hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
int			size_x, size_y, nb_steps;
double		dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
bool		file_export;
std::string export_path;

int main(int argc, char** argv)
{
	parse_args(argc, argv);
	printf("Command line options parsed\n");

	alloc();
	printf("Memory allocated\n");

	gauss_init();
	printf("State initialised\n");

	forward();
	printf("State computed\n");

	dealloc();
	printf("Memory freed\n");

	return EXIT_SUCCESS;
}
