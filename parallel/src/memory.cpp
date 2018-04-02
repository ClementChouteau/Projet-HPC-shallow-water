#include "shalw.h"
#include <stdlib.h>

int off = 0;

void alloc(void)
{
	if (block) {
		buffer_size = size_x * size_y;
	}
	else {
		off = (1 - id*size_y/p)*size_x;
		buffer_size = size_x * (size_y/p + 2);
	}

	hFil = (double*)calloc(2 * buffer_size, sizeof(double)) + off;
	uFil = (double*)calloc(2 * buffer_size, sizeof(double)) + off;
	vFil = (double*)calloc(2 * buffer_size, sizeof(double)) + off;
	hPhy = (double*)calloc(2 * buffer_size, sizeof(double)) + off;
	uPhy = (double*)calloc(2 * buffer_size, sizeof(double)) + off;
	vPhy = (double*)calloc(2 * buffer_size, sizeof(double)) + off;
}

void dealloc(void)
{
	free(hFil - off);
	free(uFil - off);
	free(vFil - off);
	free(hPhy - off);
	free(uPhy - off);
	free(vPhy - off);
}
