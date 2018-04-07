#include <mpi.h>
#include <stdio.h>

#include "shalw.h"


static MPI_File		fh;
static MPI_Request  request;
static MPI_Datatype filetype;

void create_file(void)
{
	request = 0;

	char fname[256];
	sprintf(fname, "%s/shalw_%dx%d_T%d.sav", export_path.c_str(), size_x,
			size_y, nb_steps);

	MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
				  MPI_INFO_NULL, &fh);

	if (block)
	{
		int gsizes[2] = {size_x, size_y};
		int lsizes[2] = {size_block_x, size_block_y};
		int starts[2] = {start_block_x, start_block_y};

		MPI_Type_create_subarray(2, gsizes, lsizes, starts, MPI_ORDER_C,
								 MPI_DOUBLE, &filetype);
		MPI_Type_commit(&filetype);

		MPI_File_set_view(fh, 0, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
	}
}

void export_step(int t)
{
	if (async)
	{
		MPI_Wait(&request, NULL);
		MPI_File_iwrite(fh, (void*)&HFIL(t, (block), 1), size_block, MPI_DOUBLE,
						&request);
	}
	else // sync
	{
		MPI_Offset disp =
			(t * size + start_block_y * size_block_x) * sizeof(double);
		MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native",
						  MPI_INFO_NULL);
		MPI_File_write(fh, (void*)&HFIL(t, (block), 1), size_block, MPI_DOUBLE,
					   MPI_STATUS_IGNORE);
	}
}

void finalize_export(void)
{
	if (request)
		MPI_Wait(&request, NULL);
	MPI_File_close(&fh);
}
