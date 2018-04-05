#include <mpi.h>
#include <stdio.h>


#include "shalw.h"

static MPI_File	fh;
static MPI_Offset  blocksize;
static MPI_Request request;

// void create_file(int bs)
// {
// 	blocksize = bs;
// 	request   = 0;

// 	char fname[256];
// 	sprintf(fname, "%s/shalw_%dx%d_T%d.sav", export_path.c_str(), size_x,
// 			size_y, nb_steps);

// 	MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
// 				  MPI_INFO_NULL, &fh);
// }

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

// void export_step_blocks(int t, int async)
// {
// 	int gsizes[2] = {size_x, size_y};
// 	int psizes[2] = {q, q};
// 	int lsizes[2] = {size_x / psizes[0], size_y / psizes[0]};
// 	int coords[2] = {id % psizes[0], id / psizes[1]};

// 	int starts[2] = {coords[0] * lsizes[0], coords[1] * lsizes[1]};
// 	if (t == 39)
// 		printf("id : %d -> (%d, %d)\n", id, starts[0], starts[1]);

// 	MPI_Datatype filetype;
// 	MPI_Type_create_subarray(2, gsizes, lsizes, starts, MPI_ORDER_C, MPI_DOUBLE,
// 							 &filetype);
// 	MPI_Type_commit(&filetype);

// 	MPI_Offset disp = t * size_x * size_y * sizeof(double);
// 	MPI_File_set_view(fh, disp, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

// 	if (request)
// 		MPI_Wait(&request, NULL);

// 	if (async)
// 		MPI_File_iwrite(fh, (void*)&HFIL(t, starts[0], starts[1]), blocksize,
// 						MPI_DOUBLE, &request);
// 	else
// 		MPI_File_write(fh, (void*)&HFIL(t, starts[0], starts[1]), blocksize,
// 					   MPI_DOUBLE, MPI_STATUS_IGNORE);
// }

// void finalize_export(void)
// {
// 	if (request)
// 		MPI_Wait(&request, NULL);
// 	MPI_File_close(&fh);
// }

static FILE* f;
void		 create_file(int bs)
{
	if (id != 0)
		return;
	char fname[256];

	sprintf(fname, "%s/shalw_%dx%d_T%d.sav", export_path.c_str(), size_x,
			size_y, nb_steps);

	f = fopen(fname, "w+b");
}

void export_step_blocks(int t, int async)
{
	if (id != 0)
		return;
	fwrite((void*)&HFIL(t, 0, 0), sizeof(double), size_x * size_y, f);
}

void finalize_export(void)
{
	if (id != 0)
		return;
	fclose(f);
}