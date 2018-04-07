#include <mpi.h>
#include <stdio.h>

#include "shalw.h"


static MPI_File		fh;
static MPI_Request  request;
static MPI_Datatype filetype, memtype;

void create_file(void)
{
	request = 0;

// 	char fname[256];
// 	sprintf(fname, "%s/shalw_%dx%d_T%d.sav", export_path.c_str(), size_x,
// 			size_y, nb_steps);

	MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY,
				  MPI_INFO_NULL, &fh);

	if (block)
	{
		int gsizes[2] = {size_y, size_x};
		int lsizes[2] = {size_block_y, size_block_x};
		int starts[2] = {start_block_y, start_block_x};

		// parce qu'on a intégré les bords dans le buffer mémoire, on ne peut
		// pas simplement considérer un sous tableau. La mémoire n'est pas
		// contigue de la même façon que le fichier. On créer donc un sous type
		// pour la mémoire : des lignes de size_block_x espacées de size_block_x
		// + 2 (bords)
		MPI_Datatype line;
		MPI_Type_contiguous(size_block_x, MPI_DOUBLE, &line);
		MPI_Type_create_resized(line, 0, (size_block_x + 2) * sizeof(double),
								&memtype);
		MPI_Type_commit(&memtype);

		// sous tableau pour écrire dans le fichier
		MPI_Type_create_subarray(2, gsizes, lsizes, starts, MPI_ORDER_C,
								 MPI_DOUBLE, &filetype);
		MPI_Type_commit(&filetype);

		MPI_File_set_view(fh, 0, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
	}
}

void export_step(int t)
{
	if (block)
	{
		if (async)
		{
			if (request)
				MPI_Wait(&request, NULL);
			MPI_File_iwrite(fh, (void*)&HFIL(t, 1, 1), size_block_y, memtype,
							&request);
		}
		else // sync
			MPI_File_write(fh, (void*)&HFIL(t, 1, 1), size_block_y, memtype,
						   MPI_STATUS_IGNORE);
	}
	else // bands
	{
		MPI_Offset disp =
			(t * size + start_block_y * size_block_x) * sizeof(double);
		MPI_File_set_view(fh, disp, MPI_DOUBLE, MPI_DOUBLE, "native",
						  MPI_INFO_NULL);
		if (async)
		{
			if (request)
				MPI_Wait(&request, NULL);
			MPI_File_iwrite(fh, (void*)&HFIL(t, 0, 1), size_block, MPI_DOUBLE,
							&request);
		}
		else // sync
			MPI_File_write(fh, (void*)&HFIL(t, 0, 1), size_block, MPI_DOUBLE,
						   MPI_STATUS_IGNORE);
	}
}

void export_step_blocks(int t, int async)
{
	if (id != 0)
		return;
	fwrite((void*)&HFIL(t, 0, 0), sizeof(double), size_x * size_y, f);
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