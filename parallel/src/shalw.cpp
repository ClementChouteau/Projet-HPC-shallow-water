#include "shalw.h"
#include "export.h"
#include "forward.h"
#include "init.h"
#include "memory.h"
#include "parse_args.hpp"
#include <mpi.h>
#include <stdlib.h>


#include <math.h>

#undef DEBUG

#ifdef DEBUG
#define PRINT(X) printf(X)
#else
#define PRINT(X)
#endif

double *	hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
int			size_x, size_y, nb_steps;
double		dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
bool		file_export;
std::string export_path;
int			p, id;
int			q;
bool		async;
bool		block;
int			buffer_size;

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	q = sqrt(p);

	parse_args(argc, argv);
	PRINT("Command line options parsed\n");

	if (block && (size_x % q != 0 || size_y % q != 0))
	{
		printf("Dimensions of image not compatible with number of block "
			   "workers\n");
		exit(1);
	}
	else if (!block && size_y % p != 0)
	{
		printf("Height of image not divisible by number of workers\n");
		exit(1);
	}

	alloc();
	PRINT("Memory allocated\n");

	gauss_init();
	PRINT("State initialised\n");

	if (async)
	{
		PRINT("Asynchonous mode\n");
		if (block)
			forward_blocks_async();
		else
			forward_bands_async();
	}
	else
	{
		PRINT("Synchonous mode\n");
		if (block)
			forward_blocks_sync();
		else
			forward_bands_sync();
	}
	PRINT("State computed\n");

	dealloc();
	PRINT("Memory freed\n");

	MPI_Finalize();

	return EXIT_SUCCESS;
}
