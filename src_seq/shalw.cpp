#include <stdlib.h>
#include <shalw.h>
#include <parse_args.hpp>
#include <memory.h>
#include <init.h>
#include <forward.h>
#include <export.h>
#include <sys/time.h>   /* chronometrage */

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

int main(int argc, char **argv) {
  double debut, fin;

  debut = my_gettimeofday();

  parse_args(argc, argv);
  printf("Command line options parsed\n");

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

  return EXIT_SUCCESS;
}
