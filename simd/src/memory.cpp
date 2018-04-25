#include <stdlib.h>
#include <cstring>
#include "shalw.h"

void alloc(void) {
  hFil = (double *) aligned_alloc(32, 2*size_x*size_y*sizeof(double)); // !size must be a multiple of alignment
  uFil = (double *) aligned_alloc(32, 2*size_x*size_y*sizeof(double));
  vFil = (double *) aligned_alloc(32, 2*size_x*size_y*sizeof(double));
  hPhy = (double *) aligned_alloc(32, 2*size_x*size_y*sizeof(double));
  uPhy = (double *) aligned_alloc(32, 2*size_x*size_y*sizeof(double));
  vPhy = (double *) aligned_alloc(32, 2*size_x*size_y*sizeof(double));

  std::memset(hFil, 0, 2*size_x*size_y*sizeof(double));
  std::memset(uFil, 0, 2*size_x*size_y*sizeof(double));
  std::memset(vFil, 0, 2*size_x*size_y*sizeof(double));
  std::memset(hPhy, 0, 2*size_x*size_y*sizeof(double));
  std::memset(uPhy, 0, 2*size_x*size_y*sizeof(double));
  std::memset(vPhy, 0, 2*size_x*size_y*sizeof(double));
}

void dealloc(void) {
  free(hFil);
  free(uFil);
  free(vFil);
  free(hPhy);
  free(uPhy);
  free(vPhy);
}
