#include "shalw.h"
#include <stdlib.h>

void alloc(void)
{
	hFil = (double*)calloc(2 * size_x * size_y, sizeof(double));
	uFil = (double*)calloc(2 * size_x * size_y, sizeof(double));
	vFil = (double*)calloc(2 * size_x * size_y, sizeof(double));
	hPhy = (double*)calloc(2 * size_x * size_y, sizeof(double));
	uPhy = (double*)calloc(2 * size_x * size_y, sizeof(double));
	vPhy = (double*)calloc(2 * size_x * size_y, sizeof(double));
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

// POUR UNE EVENTUELLE VERSION AVEC ALLOCATION MEMOIRE DES BLOCS
// POSSIBLE ET RAISONNABLE APRES MPI I/O

void alloc_blocks(void)
{
	// allocation des blocs
	// comportant des bords

	hFil = (double*)calloc(2 * (size_x / q + 2) * (size_y / q + 2),
						   sizeof(double));
	uFil = (double*)calloc(2 * (size_x / q + 2) * (size_y / q + 2),
						   sizeof(double));
	vFil = (double*)calloc(2 * (size_x / q + 2) * (size_y / q + 2),
						   sizeof(double));
	hPhy = (double*)calloc(2 * (size_x / q + 2) * (size_y / q + 2),
						   sizeof(double));
	uPhy = (double*)calloc(2 * (size_x / q + 2) * (size_y / q + 2),
						   sizeof(double));
	vPhy = (double*)calloc(2 * (size_x / q + 2) * (size_y / q + 2),
						   sizeof(double));

	// On d�cale pour permettre d'acc�der aux bords
	// https://stackoverflow.com/questions/3473675/are-negative-array-indexes-allowed-in-c
	hFil += size_x / q + 2 + 1;
	uFil += size_x / q + 2 + 1;
	vFil += size_x / q + 2 + 1;
	hPhy += size_x / q + 2 + 1;
	uPhy += size_x / q + 2 + 1;
	vPhy += size_x / q + 2 + 1;
}

void dealloc_blocks(void)
{
	hFil -= size_x / q + 2 + 1;
	uFil -= size_x / q + 2 + 1;
	vFil -= size_x / q + 2 + 1;
	hPhy -= size_x / q + 2 + 1;
	uPhy -= size_x / q + 2 + 1;
	vPhy -= size_x / q + 2 + 1;

	free(hFil);
	free(uFil);
	free(vFil);
	free(hPhy);
	free(uPhy);
	free(vPhy);
}
