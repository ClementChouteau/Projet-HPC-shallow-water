#include <mpi.h>
#include <stdio.h>


#include "shalw.h"

static MPI_File   fh;
static MPI_Offset blocksize;
static MPI_Offset disp;

void create_file(int bs)
{
	blocksize = bs;

	char fname[256];
	sprintf(fname, "%s/shalw_%dx%d_T%d.sav", export_path.c_str(), size_x,
			size_y, nb_steps);

	MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
				  MPI_INFO_NULL, &fh);

	// Very time-consuming operation : it's a collective method
	// MPI_File_preallocate(fh, (size_x * size_y * nb_steps) *
	// sizeof(MPI_DOUBLE));
}

void export_step_bands_sync(int t)
{
	disp = id * blocksize * sizeof(double);

	// Same displacement is applied on file and memory
	MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native",
					  MPI_INFO_NULL);

	MPI_File_write(fh, (void*)&HFIL(t, 0, 0), blocksize, MPI_DOUBLE,
				   MPI_STATUS_IGNORE);
}

void export_step_blocks_sync(int t) {}

void finalize_export(void)
{
	MPI_File_close(&fh);
}
