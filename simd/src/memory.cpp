#include <stdlib.h>
#include <cstring>
#include "shalw.h"

static double* _huvFilPhy;

void alloc(void) {
	_huv_ptr = (double *) aligned_alloc(32, 6*2*size_x*size_y*sizeof(double));

	hFil = _huvFilPhy+0*2*size_x*size_y;
	uFil = _huvFilPhy+1*2*size_x*size_y;
	vFil = _huvFilPhy+2*2*size_x*size_y;
	hPhy = _huvFilPhy+3*2*size_x*size_y;
	uPhy = _huvFilPhy+4*2*size_x*size_y;
	vPhy = _huvFilPhy+5*2*size_x*size_y;

	std::memset(hFil, 0, 6*2*size_x*size_y*sizeof(double));
}

void dealloc(void) {
	free(_huvFilPhy);
}
