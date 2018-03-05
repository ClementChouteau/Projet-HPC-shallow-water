#include <stdlib.h>
#include <shalw.h>
#include <parse_args.hpp>
#include <memory.h>
#include <init.h>
#include <forward.h>
#include <export.h>
#include <sys/time.h>   /* chronometrage */
#include <mpi.h>

double my_gettimeofday(){
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
int size_x, size_y, nb_steps;
double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
bool file_export;
std::string export_path;
int p, id;

int main(int argc, char **argv) {
	double debut, fin;

	debut = my_gettimeofday();

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	parse_args(argc, argv);
	printf("Command line options parsed\n");

	if (size_y % p != 0) {
		printf("Height of image not divisible by number of workers\n");
		exit(1);
	}

	// 1 seul export fait par root
	if (p > 1 && id != 0)
		file_export = false;

	alloc();
	printf("Memory allocated\n");

	gauss_init();
	printf("State initialised\n");

	forward();
	printf("State computed\n");

	fin = my_gettimeofday();
	printf("Temps total de calcul : %g seconde(s) \n", fin - debut);

	dealloc();
	printf("Memory freed\n");

	MPI_Finalize();

	return EXIT_SUCCESS;
}
