#include "shalw.h"
#include <math.h>

void gauss_init(void)
{
	const double gmx = size_x * dx / 2;
	const double gmy = size_y * dy / 2;
	const double gsx = 25000;
	const double gsy = 25000;

	int start_x = start_block_x;
	int start_y = start_block_y - 1; // one extra line en top
	int end_x   = end_block_x;
	int end_y   = end_block_y + 1; // one extra line en bottom

	if (block)
	{
		start_x -= 1; // one extra column on left
		end_x += 1;   // one extra column on right
	}

	for (int y = std::max(0, start_y); y < std::min(size_y, end_y); y++)
	for (int x = std::max(0, start_x); x < std::min(size_x, end_x); x++)
	{
		HFIL(0, x - start_x, y - start_y) =
			height * (exp(-pow((x * dx - gmx) / gsx, 2) / 2.)) *
			(exp(-pow((y * dy - gmy) / gsy, 2) / 2.));
	}
}
