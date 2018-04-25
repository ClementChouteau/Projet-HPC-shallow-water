#include <stdlib.h>
#include "shalw.h"
#include "parse_args.hpp"
#include "memory.h"
#include "init.h"
#include "forward.h"
#include "export.h"

double *__restrict__ hFil, *__restrict__ uFil, *__restrict__ vFil, *__restrict__ hPhy, *__restrict__ uPhy, *__restrict__ vPhy;
int size_x, size_y, nb_steps;
double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
bool file_export;
std::string export_path;

int main(int argc, char **argv) {
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
