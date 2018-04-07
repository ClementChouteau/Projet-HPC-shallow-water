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

	printf("id %d init from (%d, %d) to (%d, %d)\n", id, start_x, start_y,
		   end_x, end_y);
	for (int x = start_x; x < end_x; x++)
	{
		for (int y = start_y; y < end_y; y++)
		{
			if (x < 0 || y < 0 || x > size_x || y > size_y)
				continue;
			HFIL(0, x - start_x, y - start_y) =
				height * (exp(-pow((x * dx - gmx) / gsx, 2) / 2.)) *
				(exp(-pow((y * dy - gmy) / gsy, 2) / 2.));
		}
	}
}