#include "shalw.h"
#include <stdlib.h>

void alloc(void)
{
	if (block) // 2 extra columns + 2 extra lines
		buffer_size = (size_block_x + 2) * (size_block_y + 2);
	else // only 2 extra lines
		buffer_size = size_block_x * (size_block_y + 2);

	hFil = (double*)calloc(2 * buffer_size, sizeof(double));
	uFil = (double*)calloc(2 * buffer_size, sizeof(double));
	vFil = (double*)calloc(2 * buffer_size, sizeof(double));
	hPhy = (double*)calloc(2 * buffer_size, sizeof(double));
	uPhy = (double*)calloc(2 * buffer_size, sizeof(double));
	vPhy = (double*)calloc(2 * buffer_size, sizeof(double));
}

void dealloc(void)
{
	free(hFil);
	free(uFil);
	free(vFil);
	free(hPhy);
	free(uPhy);
	free(vPhy);
}
