#include <mpi.h>
#include <stdio.h>


#include "shalw.h"

static MPI_File	fh;
static MPI_Offset  blocksize;
static MPI_Request request;

void create_file(int bs)
{
	blocksize = bs;
	request   = 0;

	char fname[256];
	sprintf(fname, "%s/shalw_%dx%d_T%d.sav", export_path.c_str(), size_x,
			size_y, nb_steps);

	MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
				  MPI_INFO_NULL, &fh);
}

void export_step_bands(int t, int async)
{
	MPI_Offset disp = (t * p + id) * blocksize * sizeof(double);

	// Same displacement is applied on file and memory
	MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native",
					  MPI_INFO_NULL);

	if (request)
		MPI_Wait(&request, NULL);

	if (async)
		MPI_File_iwrite(fh, (void*)&HFIL(t, 0, id * (size_y / p)), blocksize,
						MPI_DOUBLE, &request);
	else
		MPI_File_write(fh, (void*)&HFIL(t, 0, id * (size_y / p)), blocksize,
					   MPI_DOUBLE, MPI_STATUS_IGNORE);
}

void export_step_blocks(int t, int async)
{
	MPI_Datatype filetype;
	int			 gsizes[2] = {size_x, size_y};
	int			 lsizes[2] = {size_x / q, size_y / q};
	int			 psizes[2] = {q, q};
	int			 coords[2] = {id % psizes[0], id / psizes[1]};
	int			 starts[2] = {coords[0] * lsizes[0], coords[1] * lsizes[1]};

	MPI_Type_create_subarray(2, gsizes, lsizes, starts, MPI_ORDER_C, MPI_DOUBLE,
							 &filetype);
	MPI_Type_commit(&filetype);

	MPI_File_set_view(fh, 0, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

	if (request)
		MPI_Wait(&request, NULL);

	if (async)
		MPI_File_iwrite(fh,
						(void*)&HFIL(t, coords[0] * (gsizes[0] / psizes[0]),
									 coords[1] * (gsizes[1] / psizes[1])),
						blocksize, MPI_DOUBLE, &request);
	else
		MPI_File_write(fh,
					   (void*)&HFIL(t, coords[0] * (gsizes[0] / psizes[0]),
									coords[1] * (gsizes[1] / psizes[1])),
					   blocksize, MPI_DOUBLE, MPI_STATUS_IGNORE);
}

void finalize_export(void)
{
	if (request)
		MPI_Wait(&request, NULL);
	MPI_File_close(&fh);
}
