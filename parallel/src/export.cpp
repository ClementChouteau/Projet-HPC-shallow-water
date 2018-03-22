#include "shalw.h"
#include <stdio.h>

FILE* create_file(void)
{
	FILE* f;
	char  fname[256];

	sprintf(fname, "%s/shalw_%dx%d_T%d.sav", export_path.c_str(), size_x,
			size_y, nb_steps);

	f = fopen(fname, "w+b");

	return f;
}

void export_step(FILE* f, int t)
{
	fwrite((void*)&HFIL(t, 0, 0), sizeof(double), size_x * size_y, f);
}

void finalize_export(FILE* f)
{
	fclose(f);
}

/*
MPI_File* create_file(void)
{
	MPI_File* f;
	char	  fname[256];

	sprintf(fname, "%s/shalw_%dx%d_T%d.sav", export_path.c_str(), size_x,
			size_y, nb_steps);

	// ATTENTION IL Y A PLUSIEURS TYPES D'IO MPI, C'EST PEUT ETRE PAS LE BON

	MPI_File_open(MPI_COMM_WORLD, "test.out", MPI_MODE_CREATE | MPI_MODE_WRONLY,
				  MPI_INFO_NULL, &f);
	MPI_File_preallocate(f, (size_x * size_y * nb_steps) * sizeof(MPI_DOUBLE));

	return f;
}

void export_step(MPI_File* f, int t)
{
	const MPI_Offset offset = (t * p + id) * ((size_x * size_y) / p);

	// CHANGER DES VALEURS SELON QUE L'ON SOIT EN MODE BLOC OU BANDES
	// EN MODE BANDES ON PEUT FAIRE UN UNIQUE WRITE
	// EN MODE BLOC ON DOIT FAIRE PLUSIEURS WRITES CAR IL Y A LES BORDS DANS LE
	// TABLEAU

	for (int j = 0; j < ; j++)
		MPI_File_write_at(
			f, offset + j * size_y,
			&HFIL(t, id_x * (size_x / q), id_y * (size_y / q) + j), size_y,
			MPI_DOUBLE, NULL);
}

void finalize_export(MPI_File* f)
{
	MPI_File_close(&h);
}
*/
