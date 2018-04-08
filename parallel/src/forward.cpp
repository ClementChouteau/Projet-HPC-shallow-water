#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "export.h"
#include "forward.h"
#include "shalw.h"


double hFil_forward(int t, int i, int j)
{
	// Phase d'initialisation du filtre
	// HPHY(t - 1, i, j) est encore nul
	if (t <= 2)
		return HPHY(t, i, j);
	return HPHY(t - 1, i, j) +
		   alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

double uFil_forward(int t, int i, int j)
{
	// Phase d'initialisation du filtre
	// UPHY(t - 1, i, j) est encore nul
	if (t <= 2)
		return UPHY(t, i, j);
	return UPHY(t - 1, i, j) +
		   alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

double vFil_forward(int t, int i, int j)
{
	// Phase d'initialisation du filtre
	// VPHY(t - 1, i, j) est encore nul
	if (t <= 2)
		return VPHY(t, i, j);
	return VPHY(t - 1, i, j) +
		   alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
}

double hPhy_forward(int t, int i, int j)
{
	double c, d;

	c = 0.;
	if (i > 0)
		c = UPHY(t - 1, i - 1, j);

	d = 0.;
	if (j < size_y - 1)
		d = VPHY(t - 1, i, j + 1);

	return HFIL(t - 1, i, j) -
		   dt * hmoy *
			   ((UPHY(t - 1, i, j) - c) / dx + (d - VPHY(t - 1, i, j)) / dy);
}

double uPhy_forward(int t, int i, int j)
{
	double b, e, f, g;

	if (i == size_x - 1)
		return 0.;

	b = 0.;
	if (i < size_x - 1)
		b = HPHY(t - 1, i + 1, j);

	e = 0.;
	if (j < size_y - 1)
		e = VPHY(t - 1, i, j + 1);

	f = 0.;
	if (i < size_x - 1)
		f = VPHY(t - 1, i + 1, j);

	g = 0.;
	if (i < size_x - 1 && j < size_y - 1)
		g = VPHY(t - 1, i + 1, j + 1);

	return UFIL(t - 1, i, j) +
		   dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) +
				 (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) -
				 (dissip * UFIL(t - 1, i, j)));
}

double vPhy_forward(int t, int i, int j)
{
	double c, d, e, f;

	if (j == 0)
		return 0.;

	c = 0.;
	if (j > 0)
		c = HPHY(t - 1, i, j - 1);

	d = 0.;
	if (i > 0 && j > 0)
		d = UPHY(t - 1, i - 1, j - 1);

	e = 0.;
	if (i > 0)
		e = UPHY(t - 1, i - 1, j);

	f = 0.;
	if (j > 0)
		f = UPHY(t - 1, i, j - 1);

	return VFIL(t - 1, i, j) +
		   dt * ((-grav / dy) * (HPHY(t - 1, i, j) - c) -
				 (pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) -
				 (dissip * VFIL(t - 1, i, j)));
}

void FORWARD(int t, int i, int j)
{
	HPHY(t, i, j) = hPhy_forward(t, i, j);
	HFIL(t, i, j) = hFil_forward(t, i, j);

	UPHY(t, i, j) = uPhy_forward(t, i, j);
	UFIL(t, i, j) = uFil_forward(t, i, j);

	VPHY(t, i, j) = vPhy_forward(t, i, j);
	VFIL(t, i, j) = vFil_forward(t, i, j);
}

// Cette fonction est longue et concentre les versions bandes/blocks
// sync/async pour factoriser le code et bien comprendre ce que ces modes
// impliquent sur une version standard. Ne pas hésiter à masquer des boucles
// ou des branchements dans un éditeur de code pour avoir un aperçu global
// de la fonction.
void forward(void)
{
	double svdt		 = 0.;
	int	t		 = 0;
	double total_msg = 0, total_calc = 0, total_export = 0;

	// Special type MPI to exchange array columns
	MPI_Datatype column;
	MPI_Type_vector(size_block_y, 1, size_block_x + 2, MPI_DOUBLE, &column);
	MPI_Type_commit(&column);

	// Requests for async mode
	MPI_Request r[8], s[8];
	for (int i = 0; i < 8; i++)
	{
		r[i] = MPI_REQUEST_NULL;
		s[i] = MPI_REQUEST_NULL;
	}

	if (file_export)
		create_file();

	// t = 0 is the initial state already in memory by gauss_init
	for (t = 1; t < nb_steps; t++) // we talk about t iterations but actually we
								   // calculate only t - 1 times
	{
		if (t == 1)
		{
			svdt = dt;
			dt   = 0;
		}
		if (t == 2)
			dt = svdt / 2.;

		if (file_export && t % step_export == 0)
		{
			clock_t start_export = clock();
			export_step(t - 1); // t - 1 is ready to export
			total_export += TIME(start_export, clock());
		}
		// recouvrement par le calcul en async

		// ECHANGE DE LIGNES
		// à l'instant t, la grille t-1 est calculée, on peut donc échanger les
		// lignes de t - 1
		clock_t start_msg = clock();
		if (t > 1)
		{
			if (async)
			{
				if (block)
				{
					// Echanges des Colonnes
					if (id_x > 0) // not first column blocks process
					{
						// Send first column HPHY to id_x - 1
						MPI_Isend(&HPHY(t - 1, 1, 1), 1, column,
								  ID(id_x - 1, id_y), 0, MPI_COMM_WORLD, s);

						// Receive first column UPHY from id_x - 1
						MPI_Irecv(&UPHY(t - 1, 0, 1), 1, column,
								  ID(id_x - 1, id_y), 0, MPI_COMM_WORLD, r);

						// Send first column VPHY to id_x - 1
						MPI_Isend(&VPHY(t - 1, 1, 1), 1, column,
								  ID(id_x - 1, id_y), 0, MPI_COMM_WORLD, s + 1);
					}
					if (id_x < p_x - 1) // not last column blocks process
					{
						// Receive last column HPHY from id_x + 1
						MPI_Irecv(&HPHY(t - 1, size_block_x + 1, 1), 1, column,
								  ID(id_x + 1, id_y), 0, MPI_COMM_WORLD, r + 1);

						// Send last column UPHY to id_x + 1
						MPI_Isend(&UPHY(t - 1, size_block_x, 1), 1, column,
								  ID(id_x + 1, id_y), 0, MPI_COMM_WORLD, s + 2);

						// Receive last column VPHY from id_x + 1
						MPI_Irecv(&VPHY(t - 1, size_block_x + 1, 1), 1, column,
								  ID(id_x + 1, id_y), 0, MPI_COMM_WORLD, r + 2);
					}

					// Echanges des lignes
					if (id_y > 0) // not first line blocks process
					{
						// Receive first line HPHY from id_y - 1
						MPI_Irecv(&HPHY(t - 1, 1, 0), size_block_x, MPI_DOUBLE,
								  ID(id_x, id_y - 1), 0, MPI_COMM_WORLD, r + 3);

						// Receive first line UPHY from id_y - 1
						MPI_Irecv(&UPHY(t - 1, 1, 0), size_block_x, MPI_DOUBLE,
								  ID(id_x, id_y - 1), 0, MPI_COMM_WORLD, r + 4);

						// Send first line VPHY to id_y - 1
						MPI_Isend(&VPHY(t - 1, 1, 1), size_block_x, MPI_DOUBLE,
								  ID(id_x, id_y - 1), 0, MPI_COMM_WORLD, s + 3);
					}
					if (id_y < p_y - 1) // not last line blocks process
					{
						// Send last line HPHY to id_y + 1
						MPI_Isend(&HPHY(t - 1, 1, size_block_y), size_block_x,
								  MPI_DOUBLE, ID(id_x, id_y + 1), 0,
								  MPI_COMM_WORLD, s + 4);

						// Send last line UPHY to id_y + 1
						MPI_Isend(&UPHY(t - 1, 1, size_block_y), size_block_x,
								  MPI_DOUBLE, ID(id_x, id_y + 1), 0,
								  MPI_COMM_WORLD, s + 5);

						// Receive last line VPHY from id_y + 1
						MPI_Irecv(&VPHY(t - 1, 1, size_block_y + 1),
								  size_block_x, MPI_DOUBLE, ID(id_x, id_y + 1),
								  0, MPI_COMM_WORLD, r + 5);
					}

					// Echanges des coins
					if (id_x > 0 &&
						id_y > 0) // not first top left block process
					{
						// Receive first top left UPHY value from (id_x - 1,
						// id_y - 1)
						MPI_Irecv(&UPHY(t - 1, 0, 0), 1, MPI_DOUBLE,
								  ID(id_x - 1, id_y - 1), 0, MPI_COMM_WORLD,
								  r + 6);

						// Send first top left VPHY value to (id_x - 1, id_y -
						// 1)
						MPI_Isend(&VPHY(t - 1, 1, 1), 1, MPI_DOUBLE,
								  ID(id_x - 1, id_y - 1), 0, MPI_COMM_WORLD,
								  s + 6);
					}
					if (id_x < p_x - 1 &&
						id_y < p_y - 1) // not last bottom right block process
					{
						// Send last bottom right UPHY value to (id_x + 1, id_y
						// + 1)
						MPI_Isend(&UPHY(t - 1, size_block_x, size_block_y), 1,
								  MPI_DOUBLE, ID(id_x + 1, id_y + 1), 0,
								  MPI_COMM_WORLD, s + 7);

						// Receive last bottom right UPHY value from (id_x + 1,
						// id_y
						// + 1)
						MPI_Irecv(
							&VPHY(t - 1, size_block_x + 1, size_block_y + 1), 1,
							MPI_DOUBLE, ID(id_x + 1, id_y + 1), 0,
							MPI_COMM_WORLD, r + 7);
					}
				}
				else // bands
				{
					// Echange id-1 <=> id
					if (id > 0)
					{
						MPI_Irecv(&HPHY(t - 1, 0, 0), size_block_x, MPI_DOUBLE,
								  id - 1, 0, MPI_COMM_WORLD, r);
						MPI_Irecv(&UPHY(t - 1, 0, 0), size_block_x, MPI_DOUBLE,
								  id - 1, 0, MPI_COMM_WORLD, r + 1);
						MPI_Isend(&VPHY(t - 1, 0, 1), size_block_x, MPI_DOUBLE,
								  id - 1, 0, MPI_COMM_WORLD, s);
					}
					// Echange id <=> id+1
					if (id < p - 1)
					{
						MPI_Isend(&HPHY(t - 1, 0, size_block_y), size_block_x,
								  MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD, s + 1);
						MPI_Isend(&UPHY(t - 1, 0, size_block_y), size_block_x,
								  MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD, s + 2);
						MPI_Irecv(&VPHY(t - 1, 0, size_block_y + 1),
								  size_block_x, MPI_DOUBLE, id + 1, 0,
								  MPI_COMM_WORLD, r + 2);
					}
				}
			}
			else // sync
			{
				if (block)
				{
					// Echanges des Colonnes
					if (id_x > 0) // not first column blocks process
					{
						// Send first column HPHY to id_x - 1
						MPI_Send(&HPHY(t - 1, 1, 1), 1, column,
								 ID(id_x - 1, id_y), 0, MPI_COMM_WORLD);

						// Receive first column UPHY from id_x - 1
						MPI_Recv(&UPHY(t - 1, 0, 1), 1, column,
								 ID(id_x - 1, id_y), 0, MPI_COMM_WORLD, NULL);

						// Send first column VPHY to id_x - 1
						MPI_Send(&VPHY(t - 1, 1, 1), 1, column,
								 ID(id_x - 1, id_y), 0, MPI_COMM_WORLD);
					}
					if (id_x < p_x - 1) // not last column blocks process
					{
						// Receive last column HPHY from id_x + 1
						MPI_Recv(&HPHY(t - 1, size_block_x + 1, 1), 1, column,
								 ID(id_x + 1, id_y), 0, MPI_COMM_WORLD, NULL);

						// Send last column UPHY to id_x + 1
						MPI_Send(&UPHY(t - 1, size_block_x, 1), 1, column,
								 ID(id_x + 1, id_y), 0, MPI_COMM_WORLD);

						// Receive last column VPHY from id_x + 1
						MPI_Recv(&VPHY(t - 1, size_block_x + 1, 1), 1, column,
								 ID(id_x + 1, id_y), 0, MPI_COMM_WORLD, NULL);
					}

					// Echanges des lignes
					if (id_y > 0) // not first line blocks process
					{
						// Receive first line HPHY from id_y - 1
						MPI_Recv(&HPHY(t - 1, 1, 0), size_block_x, MPI_DOUBLE,
								 ID(id_x, id_y - 1), 0, MPI_COMM_WORLD, NULL);

						// Receive first line UPHY from id_y - 1
						MPI_Recv(&UPHY(t - 1, 1, 0), size_block_x, MPI_DOUBLE,
								 ID(id_x, id_y - 1), 0, MPI_COMM_WORLD, NULL);

						// Send first line VPHY to id_y - 1
						MPI_Send(&VPHY(t - 1, 1, 1), size_block_x, MPI_DOUBLE,
								 ID(id_x, id_y - 1), 0, MPI_COMM_WORLD);
					}
					if (id_y < p_y - 1) // not last line blocks process
					{
						// Send last line HPHY to id_y + 1
						MPI_Send(&HPHY(t - 1, 1, size_block_y), size_block_x,
								 MPI_DOUBLE, ID(id_x, id_y + 1), 0,
								 MPI_COMM_WORLD);

						// Send last line UPHY to id_y + 1
						MPI_Send(&UPHY(t - 1, 1, size_block_y), size_block_x,
								 MPI_DOUBLE, ID(id_x, id_y + 1), 0,
								 MPI_COMM_WORLD);

						// Receive last line VPHY from id_y + 1
						MPI_Recv(&VPHY(t - 1, 1, size_block_y + 1),
								 size_block_x, MPI_DOUBLE, ID(id_x, id_y + 1),
								 0, MPI_COMM_WORLD, NULL);
					}

					// Echanges des coins
					if (id_x > 0 &&
						id_y > 0) // not first top left block process
					{
						// Receive first top left UPHY value from (id_x - 1,
						// id_y - 1)
						MPI_Recv(&UPHY(t - 1, 0, 0), 1, MPI_DOUBLE,
								 ID(id_x - 1, id_y - 1), 0, MPI_COMM_WORLD,
								 NULL);

						// Send first top left VPHY value to (id_x-1, id_y-1)
						MPI_Send(&VPHY(t - 1, 1, 1), 1, MPI_DOUBLE,
								 ID(id_x - 1, id_y - 1), 0, MPI_COMM_WORLD);
					}
					if (id_x < p_x - 1 &&
						id_y < p_y - 1) // not last bottom right block process
					{
						// Send last bottom right UPHY value to (id_x + 1, id_y
						// + 1)
						MPI_Send(&UPHY(t - 1, size_block_x, size_block_y), 1,
								 MPI_DOUBLE, ID(id_x + 1, id_y + 1), 0,
								 MPI_COMM_WORLD);

						// Receive last bottom right UPHY value from (id_x + 1,
						// id_y + 1)
						MPI_Recv(
							&VPHY(t - 1, size_block_x + 1, size_block_y + 1), 1,
							MPI_DOUBLE, ID(id_x + 1, id_y + 1), 0,
							MPI_COMM_WORLD, NULL);
					}
				}
				else // bands
				{
					// Echange id-1 <=> id
					if (id > 0)
					{
						MPI_Recv(&HPHY(t - 1, 0, 0), size_block_x, MPI_DOUBLE,
								 id - 1, 0, MPI_COMM_WORLD, NULL);
						MPI_Recv(&UPHY(t - 1, 0, 0), size_block_x, MPI_DOUBLE,
								 id - 1, 0, MPI_COMM_WORLD, NULL);
						MPI_Send(&VPHY(t - 1, 0, 1), size_block_x, MPI_DOUBLE,
								 id - 1, 0, MPI_COMM_WORLD);
					}

					// Echange id <=> id+1
					if (id < p - 1)
					{
						MPI_Send(&HPHY(t - 1, 0, size_block_y), size_block_x,
								 MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD);
						MPI_Send(&UPHY(t - 1, 0, size_block_y), size_block_x,
								 MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD);
						MPI_Recv(&VPHY(t - 1, 0, size_block_y + 1),
								 size_block_x, MPI_DOUBLE, id + 1, 0,
								 MPI_COMM_WORLD, NULL);
					}
				}
			}
		}
		total_msg += TIME(start_msg, clock());

		// CALCULATIONS PREPARATION
		int start_x = 0;
		int start_y = 1; // skip first extra line
		int end_x   = size_block_x;
		int end_y   = size_block_y + 1; // skip last extra line

		if (block)
		{
			start_x += 1; // skip first extra column
			end_x += 1;   // skip last extra column
		}

		if (async) // we have to reduce block by 1 line for each sides
		{
			start_y += 1; // no calculations for first line
			if (block)
			{
				start_x += 1; // no calculations for first column
				end_x -= 1;   // no calculations for last colomn
			}
			end_y -= 1; // no calculations for last line
		}

		// HERE ARE MOST CALCULATIONS for t
		// Peut facilement être parallélisé avec OpenMP
		// if async mode, messages are exchanged at the same time
		clock_t start_calc = clock();
		for (int y = start_y; y < end_y; y++)
			for (int x = start_x; x < end_x; x++)
				FORWARD(t, x, y);
		total_calc += TIME(start_calc, clock());

		if (async) // Vérifier échange des bords t-1 avant de finir les calculs
		{
			// We need these receptions before finish calculations
			// Should already be finished
			clock_t start_msg = clock();
			MPI_Waitall(8, r, MPI_STATUSES_IGNORE); // Attente réception bords
			total_msg += TIME(start_msg, clock());

			// On fini les calculs des bords
			start_x -= 1;
			end_x += 1;
			clock_t start_calc = clock();
			for (int x = start_x; x < end_x; x++)
			{
				FORWARD(t, x, 1);			 // first calculation line
				FORWARD(t, x, size_block_y); // last calculation line
			}
			if (block)
			{
				for (int y = 1; y < size_block_y; y++)
				{
					FORWARD(t, 1, y);			 // first calculation column
					FORWARD(t, size_block_x, y); // last calculation column
				}
			}
			total_calc += TIME(start_calc, clock());

			// No need to wait before for these before calculations
			start_msg = clock();
			MPI_Waitall(8, s, MPI_STATUSES_IGNORE); // Attente envoi bords
			total_msg += TIME(start_msg, clock());
			// All messages have been exchanged : we can start new ones
		}

		if (t == 2)
			dt = svdt;
	}

	if (file_export)
	{
		clock_t start_export = clock();
		export_step(t - 1); // final iteration ready to export
		finalize_export();
		total_export += TIME(start_export, clock());
	}

	ID0_(printf("	Message exchange : %.2f\n", total_msg))
	ID0_(printf("	Calculations : %.2f\n", total_calc))
	if (file_export)
	{
		ID0_(printf("	Export : %.2f\n", total_export))
	}

	MPI_Barrier(MPI_COMM_WORLD);
}
