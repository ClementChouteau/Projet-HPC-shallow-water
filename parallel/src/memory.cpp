#include "shalw.h"
#include <cstdlib>

void alloc(void)
{
	if (block) // 2 extra columns + 2 extra lines
		buffer_size = (size_block_x + 2) * (size_block_y + 2);
	else // only 2 extra lines
		buffer_size = size_block_x * (size_block_y + 2);

	hFil = (double*)std::calloc(2 * buffer_size, sizeof(double));
	uFil = (double*)std::calloc(2 * buffer_size, sizeof(double));
	vFil = (double*)std::calloc(2 * buffer_size, sizeof(double));
	hPhy = (double*)std::calloc(2 * buffer_size, sizeof(double));
	uPhy = (double*)std::calloc(2 * buffer_size, sizeof(double));
	vPhy = (double*)std::calloc(2 * buffer_size, sizeof(double));
}

void dealloc(void)
{
	std::free(hFil);
	std::free(uFil);
	std::free(vFil);
	std::free(hPhy);
	std::free(uPhy);
	std::free(vPhy);
}
