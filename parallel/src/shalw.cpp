#include "shalw.h"
#include "export.h"
#include "forward.h"
#include "init.h"
#include "memory.h"
#include "parse_args.hpp"
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include <math.h>

double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
int		size_x, size_y, size, nb_steps, size_block_x, size_block_y, size_block;
int		start_block_x, start_block_y, end_block_x, end_block_y;
double  dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
bool	file_export;
std::string export_path;
int			p, id, id_x, id_y, p_x, p_y;
bool		async;
bool		block;
int			buffer_size;

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	parse_args(argc, argv);
	if (id == 0)
		printf("Command line options parsed\n");

	if (block)
	{
		p_x = p_y = sqrt(p);
		id_x	  = id % p_x;
		id_y	  = id / p_y;
	}
	else // bands
	{
		p_x  = 1;
		p_y  = p;
		id_x = 0;
		id_y = id;
	}
	size		  = size_x * size_y;
	size_block_x  = size_x / p_x;
	size_block_y  = size_y / p_y;
	size_block	= size_block_x * size_block_y;
	start_block_x = id_x * size_block_x;
	start_block_y = id_y * size_block_y;

	if (block)
	{
		end_block_x = (id_x + 1) % (p_x + 1) * size_block_x;
		end_block_y = (id_y + 1) % (p_y + 1) * size_block_y;
	}
	else // bands
	{
		end_block_x = size_block_x;
		end_block_y = (id_y + 1) * size_block_y;
	}
	// Theses variables doesn't take into account the extra columns / lines
	// needed for calculations.

	if (block && (size_x % p_x != 0 || size_y % p_y != 0))
	{
		if (id == 0)
			printf("Dimensions of image not compatible with number of block "
				   "workers\n");
		exit(1);
	}
	else if (!block && size_y % p != 0)
	{
		if (id == 0)
			printf("Height of image not divisible by number of workers\n");
		exit(1);
	}

	alloc();
	if (id == 0)
		printf("Memory allocated\n");

	gauss_init();

	if (id == 0)
	{
		printf("State initialised\n");
		printf("%s mode\n", (block) ? "Blocks" : "Bands");
		printf("%s mode\n", (async) ? "Asynchonous" : "Synchronus");
	}

	forward();

	if (id == 0)
		printf("State computed\n");

	dealloc();
	if (id == 0)
		printf("Memory freed\n");

	MPI_Finalize();

	return EXIT_SUCCESS;
}
