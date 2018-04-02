#include "shalw.h"
#include <math.h>

void gauss_init(void)
{
	double gmx, gmy, gsx, gsy;

	gmx = size_x * dx / 2;
	gmy = size_y * dy / 2;
	gsx = 25000;
	gsy = 25000;

	for (int i = 0; i < size_x; i++)
	{
		if (block) {
			for (int j = 0; j < size_y; j++)
			{
				HFIL(0, i, j) = height * (exp(-pow((i * dx - gmx) / gsx, 2) / 2.)) *
								(exp(-pow((j * dy - gmy) / gsy, 2) / 2.));
			}
		}
		else {
			for (int j = id * (size_y / p)-1; j < (id + 1) * (size_y / p)+1; j++)
			{
				HFIL(0, i, j) = height * (exp(-pow((i * dx - gmx) / gsx, 2) / 2.)) *
								(exp(-pow((j * dy - gmy) / gsy, 2) / 2.));
			}
		}
	}
}
