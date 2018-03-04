#include <stdlib.h>
#include <shalw.h>

void alloc(void) {
  hFil = (double *) calloc(2*size_x*size_y, sizeof(double));
  uFil = (double *) calloc(2*size_x*size_y, sizeof(double));
  vFil = (double *) calloc(2*size_x*size_y, sizeof(double));
  hPhy = (double *) calloc(2*size_x*size_y, sizeof(double));
  uPhy = (double *) calloc(2*size_x*size_y, sizeof(double));
  vPhy = (double *) calloc(2*size_x*size_y, sizeof(double));
}

void dealloc(void) {
  free(hFil);
  free(uFil);
  free(vFil);
  free(hPhy);
  free(uPhy);
  free(vPhy);
}
