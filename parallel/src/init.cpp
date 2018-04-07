#include "shalw.h"
#include <math.h>

void gauss_init(void)
{
	double gmx = size_x * dx / 2;
	double gmy = size_y * dy / 2;
	double gsx = 25000;
	double gsy = 25000;

	int start_x = start_block_x;
	int start_y = start_block_y - 1; // one extra line en top
	int end_x   = end_block_x;
	int end_y   = end_block_y + 1; // one extra line en bottom

	if (block)
	{
		start_x -= 1; // one extra column on left
		end_x += 1;   // one extra column on right
	}

	for (int x = start_x; x < end_x; x++)
	{
		for (int y = start_y; y < end_y; y++)
		{
			HFIL(0, x - start_x, y - start_y) =
				height * (exp(-pow((x * dx - gmx) / gsx, 2) / 2.)) *
				(exp(-pow((y * dy - gmy) / gsy, 2) / 2.));
		}
	}
}