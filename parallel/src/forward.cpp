#include <stdio.h>
#include <math.h>
#include "shalw.h"
#include "export.h"
#include <stdlib.h>

#include <mpi.h>

#define U_TAG (1)
#define V_TAG (2)
#define H_TAG (3)

#define U_RECEIVED (1)
#define V_RECEIVED (2)
#define H_RECEIVED (4)

#define U_COMPUTED (1)
#define V_COMPUTED (2)
#define H_COMPUTED (4)

#define LINE(line, k) (line+k*(size_y/q))

double hFil_forward(int t, int i, int j) {
	//Phase d'initialisation du filtre
	//HPHY(t - 1, i, j) est encore nul
	if (t <= 2)
		return HPHY(t, i, j);
	return HPHY(t - 1, i, j) +
			alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

double uFil_forward(int t, int i, int j) {
	//Phase d'initialisation du filtre
	//UPHY(t - 1, i, j) est encore nul
	if (t <= 2)
		return UPHY(t, i, j);
	return UPHY(t - 1, i, j) +
			alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

double vFil_forward(int t, int i, int j) {
	//Phase d'initialisation du filtre
	//VPHY(t - 1, i, j) est encore nul
	if (t <= 2)
		return VPHY(t, i, j);
	return VPHY(t - 1, i, j) +
			alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
}

double hPhy_forward(int t, int i, int j) {
	double c, d;

	c = 0.;
	if (i > 0)
		c = UPHY(t - 1, i - 1, j);

	d = 0.;
	if (j < size_y - 1)
		d = VPHY(t - 1, i, j + 1);

	return HFIL(t - 1, i, j) -
			dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx +
						 (d - VPHY(t - 1, i, j)) / dy);
}

double uPhy_forward(int t, int i, int j) {
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

double vPhy_forward(int t, int i, int j) {
	double c, d, e, f;

	if (j == 0)
		return 0.;

	c = 0.;
	if (j > 0)
		c = HPHY(t - 1, i, j - 1);

	d = 0.;
	if (i > 0 && j > 0)
		d = UPHY(t - 1, i -1, j -1);

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

void FORWARD(int t, int i, int j) {
	HPHY(t, i, j) = hPhy_forward(t, i, j);
	HFIL(t, i, j) = hFil_forward(t, i, j);

	UPHY(t, i, j) = uPhy_forward(t, i, j);
	UFIL(t, i, j) = uFil_forward(t, i, j);

	VPHY(t, i, j) = vPhy_forward(t, i, j);
	VFIL(t, i, j) = vFil_forward(t, i, j);
}

void forward_bands_sync(void) {
	FILE *file = NULL;
	double svdt = 0.;
	int t = 0;

	if (file_export) {
		file = create_file();
		export_step(file, t);
	}

	for (t = 1; t < nb_steps; t++) {
		if (t == 1) {
			svdt = dt;
			dt = 0;
		}
		if (t == 2){
			dt = svdt / 2.;
		}

		// Travail par bandes, 1 bande par processus
		for (int j = id*(size_y/p); j < (id+1)*(size_y/p); j++)
		for (int i = 0; i < size_x; i++)
			FORWARD(t, i, j);

		// Echange id-1 <=> id
		if (id != 0) {
			MPI_Recv(&HPHY(t, 0, id*(size_y/p)-1), size_x, MPI_DOUBLE, id-1, 0, MPI_COMM_WORLD, NULL);
			MPI_Recv(&UPHY(t, 0, id*(size_y/p)-1), size_x, MPI_DOUBLE, id-1, 0, MPI_COMM_WORLD, NULL);
			MPI_Send(&VPHY(t, 0, id*(size_y/p)), size_x, MPI_DOUBLE, id-1, 0, MPI_COMM_WORLD);
		}

		// Echange id <=> id+1
		if (id != p-1) {
			MPI_Send(&HPHY(t, 0, (id+1)*(size_y/p)-1), size_x, MPI_DOUBLE, id+1, 0, MPI_COMM_WORLD);
			MPI_Send(&UPHY(t, 0, (id+1)*(size_y/p)-1), size_x, MPI_DOUBLE, id+1, 0, MPI_COMM_WORLD);
			MPI_Recv(&VPHY(t, 0, (id+1)*(size_y/p)), size_x, MPI_DOUBLE, id+1, 0, MPI_COMM_WORLD, NULL);
		}

		if (file_export) {
			export_step(file, t);
		}

		if (t == 2) {
			dt = svdt;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (file_export) {
		finalize_export(file);
	}
}




void forward_bands_async(void) {
	FILE *file = NULL;
	double svdt = 0.;
	int t = 0;

	if (file_export) {
		file = create_file();
		export_step(file, t);
	}

	MPI_Request r[3] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
	int received = U_RECEIVED|V_RECEIVED|H_RECEIVED;
	for (t = 1; t < nb_steps; t++) {
		//printf("id: %d, t: %d\n", id, t);
		if (t == 1) {
			svdt = dt;
			dt = 0;
		}
		if (t == 2){
			dt = svdt / 2.;
		}

		// Les processeurs de rang 0 et p-1 n'ont pas besoin de tout
		if (id == 0)
			received |= U_RECEIVED|H_RECEIVED;
		if (id == p-1)
			received |= V_RECEIVED;

		// Calculs Blocs(t) -> Blocs(t+1)
		MPI_Request s[3] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
		int computed = 0;
		while (computed != (U_COMPUTED|V_COMPUTED|H_COMPUTED)) {
			// Attendre d'avoir une r√©ception
			int indx;
			MPI_Waitany(3, r, &indx, MPI_STATUS_IGNORE);

			//  mettre le flag recv_flags correspondant
			switch (indx) {
			case 0: received |= U_RECEIVED; r[0] = MPI_REQUEST_NULL; break;
			case 1: received |= V_RECEIVED; r[1] = MPI_REQUEST_NULL; break;
			case 2: received |= H_RECEIVED; r[2] = MPI_REQUEST_NULL; break;
			}

			// On regarde ce que l'on peut calculer
			if (received & (U_RECEIVED|V_RECEIVED) && (computed & H_COMPUTED) == 0) {
					// Travail par bandes, 1 bande par processus
					for (int j = id*(size_y/p); j < (id+1)*(size_y/p); j++) {
					for (int i = 0; i < size_x; i++) {
						HPHY(t, i, j) = hPhy_forward(t, i, j);}}

					// envoyer ma premi√®re ligne de HPHY
					if (id != p-1)
						MPI_Isend(&HPHY(t, 0, (id+1)*(size_y/p)-1), size_x, MPI_DOUBLE, id+1, H_TAG, MPI_COMM_WORLD, &s[2]);

					// recevoir leur derni√®re ligne de HPHY
					if (id != 0)
						MPI_Irecv(&HPHY(t, 0, id*(size_y/p)-1), size_x, MPI_DOUBLE, id-1, H_TAG, MPI_COMM_WORLD, &r[2]);

					for (int j = id*(size_y/p); j < (id+1)*(size_y/p); j++) {
					for (int i = 0; i < size_x; i++) {
						HFIL(t, i, j) = hFil_forward(t, i, j);}}

					computed |= H_COMPUTED;
			}

			if (received & (V_RECEIVED|H_RECEIVED) && (computed & U_COMPUTED) == 0) {
					// Travail par bandes, 1 bande par processus
					for (int j = id*(size_y/p); j < (id+1)*(size_y/p); j++) {
					for (int i = 0; i < size_x; i++) {
						UPHY(t, i, j) = uPhy_forward(t, i, j);}}

					// envoyer ma premi√®re ligne de UPHY
					if (id != p-1)
						MPI_Isend(&UPHY(t, 0, (id+1)*(size_y/p)-1), size_x, MPI_DOUBLE, id+1, U_TAG, MPI_COMM_WORLD, &s[0]);

					// recevoir leur derni√®re ligne de UPHY
					if (id != 0)
						MPI_Irecv(&UPHY(t, 0, id*(size_y/p)-1), size_x, MPI_DOUBLE, id-1, U_TAG, MPI_COMM_WORLD, &r[0]);

					for (int j = id*(size_y/p); j < (id+1)*(size_y/p); j++) {
					for (int i = 0; i < size_x; i++) {
						UFIL(t, i, j) = uFil_forward(t, i, j);}}

					computed |= U_COMPUTED;
			}

			if (received & (U_RECEIVED|H_RECEIVED) && (computed & V_COMPUTED) == 0) {
					// Travail par bandes, 1 bande par processus
					for (int j = id*(size_y/p); j < (id+1)*(size_y/p); j++) {
					for (int i = 0; i < size_x; i++) {
						VPHY(t, i, j) = vPhy_forward(t, i, j);}}

					// envoyer ma premi√®re ligne de VPHY
					if (id != 0)
						MPI_Isend(&VPHY(t, 0, id*(size_y/p)), size_x, MPI_DOUBLE, id-1, V_TAG, MPI_COMM_WORLD, &s[1]);

					// recevoir leur derni√®re ligne de VPHY
					if (id != p-1)
						MPI_Irecv(&VPHY(t, 0, (id+1)*(size_y/p)), size_x, MPI_DOUBLE, id+1, V_TAG, MPI_COMM_WORLD, &r[1]);

					for (int j = id*(size_y/p); j < (id+1)*(size_y/p); j++) {
					for (int i = 0; i < size_x; i++) {
						VFIL(t, i, j) = vFil_forward(t, i, j);}}

					computed |= V_COMPUTED;
			}
		}

		received = 0;

		// Attente fin des envois
		MPI_Waitall(3, s, MPI_STATUSES_IGNORE);

		if (file_export) {
			export_step(file, t);
		}

		if (t == 2) {
			dt = svdt;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (file_export) {
		finalize_export(file);
	}
}


int id_from_xy(int x, int y) {
	return y*q + x;
}

void forward_blocks_sync(void) {
	FILE *file = NULL;
	double svdt = 0.;
	int t = 0;

	if (file_export) {
		file = create_file();
		export_step(file, t);
	}

	double* line = (double*) malloc((size_y/q)*sizeof(double));

	for (t = 1; t < nb_steps; t++) {
		if (t == 1) {
			svdt = dt;
			dt = 0;
		}
		if (t == 2){
			dt = svdt / 2.;
		}

		// Travail par blocks, 1 block par processus

		const int id_x = id % q;
		const int id_y = id / q;

		for (int j = id_y*(size_y/q); j < (id_y+1)*(size_y/q); j++)
		for (int i = id_x*(size_x/q); i < (id_x+1)*(size_x/q); i++)
			FORWARD(t, i, j);

		// Echanges des lignes
		if (id_y != q-1) {
			MPI_Send(&HPHY(t, id_x*(size_x/q), (id_y+1)*(size_y/q)-1), size_x/q, MPI_DOUBLE, id_from_xy(id_x, id_y+1), 0, MPI_COMM_WORLD);
			MPI_Send(&UPHY(t, id_x*(size_x/q), (id_y+1)*(size_y/q)-1), size_x/q, MPI_DOUBLE, id_from_xy(id_x, id_y+1), 0, MPI_COMM_WORLD);

			MPI_Recv(&VPHY(t, id_x*(size_x/q), (id_y+1)*(size_y/q)), size_x/q, MPI_DOUBLE, id_from_xy(id_x, id_y+1), 0, MPI_COMM_WORLD, NULL);
		}
		if (id_y != 0) {
			MPI_Recv(&HPHY(t, id_x*(size_x/q), id_y*(size_y/q)-1), size_x/q, MPI_DOUBLE, id_from_xy(id_x, id_y-1), 0, MPI_COMM_WORLD, NULL);
			MPI_Recv(&UPHY(t, id_x*(size_x/q), id_y*(size_y/q)-1), size_x/q, MPI_DOUBLE, id_from_xy(id_x, id_y-1), 0, MPI_COMM_WORLD, NULL);

			MPI_Send(&VPHY(t, id_x*(size_x/q), id_y*(size_y/q)), size_x/q, MPI_DOUBLE, id_from_xy(id_x, id_y-1), 0, MPI_COMM_WORLD);
		}

		// Echanges des Colonnes
		if (id_x != 0) {
			for (int j = 0; j < size_y/q; j++)
				line[j] = HPHY(t, id_x*(size_x/q), id_y*(size_y/q)+j);
			MPI_Send(line, size_y/q, MPI_DOUBLE, id_from_xy(id_x-1, id_y), 0, MPI_COMM_WORLD);

			for (int j = 0; j < size_y/q; j++)
				line[j] = VPHY(t, id_x*(size_x/q), id_y*(size_y/q)+j);
			MPI_Send(line, size_y/q, MPI_DOUBLE, id_from_xy(id_x-1, id_y), 0, MPI_COMM_WORLD);

			MPI_Recv(line, size_y/q, MPI_DOUBLE, id_from_xy(id_x-1, id_y), 0, MPI_COMM_WORLD, NULL);
			for (int j = 0; j < size_y/q; j++)
				UPHY(t, id_x*(size_x/q)-1, id_y*(size_y/q)+j) = line[j];
		}
		if (id_x != q-1) {
			MPI_Recv(line, size_y/q, MPI_DOUBLE, id_from_xy(id_x+1, id_y), 0, MPI_COMM_WORLD, NULL);
			for (int j = 0; j < size_y/q; j++)
				HPHY(t, (id_x+1)*(size_x/q), id_y*(size_y/q)+j) = line[j];

			MPI_Recv(line, size_y/q, MPI_DOUBLE, id_from_xy(id_x+1, id_y), 0, MPI_COMM_WORLD, NULL);
			for (int j = 0; j < size_y/q; j++)
				VPHY(t, (id_x+1)*(size_x/q), id_y*(size_y/q)+j) = line[j];

			for (int j = 0; j < size_y/q; j++)
				line[j] = UPHY(t, (id_x+1)*(size_x/q)-1, id_y*(size_y/q)+j);
			MPI_Send(line, size_y/q, MPI_DOUBLE, id_from_xy(id_x+1, id_y), 0, MPI_COMM_WORLD);
		}

		// Echanges des coins
		if (id_y != 0 && id_x != 0) {
			MPI_Recv(&UPHY(t, id_x*(size_x/q)-1, id_y*(size_y/q))-1, 1, MPI_DOUBLE, id_from_xy(id_x-1, id_y-1), 0, MPI_COMM_WORLD, NULL);
			MPI_Send(&VPHY(t, id_x*(size_x/q), id_y*(size_y/q)), 1, MPI_DOUBLE, id_from_xy(id_x-1, id_y-1), 0, MPI_COMM_WORLD);
		}

		if (id_y != q-1 && id_x != q-1) {
			MPI_Send(&UPHY(t, (id_x+1)*(size_x/q)-1, (id_y+1)*(size_y/q))-1, 1, MPI_DOUBLE, id_from_xy(id_x+1, id_y+1), 0, MPI_COMM_WORLD);
			MPI_Recv(&VPHY(t, (id_x+1)*(size_x/q), (id_y+1)*(size_y/q)), 1, MPI_DOUBLE, id_from_xy(id_x+1, id_y+1), 0, MPI_COMM_WORLD, NULL);
		}


		if (file_export) {
			export_step(file, t);
		}

		if (t == 2) {
			dt = svdt;
		}
	}

	free(line);

	if (file_export) {
		finalize_export(file);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void forward_blocks_async(void) {
	FILE *file = NULL;
	double svdt = 0.;
	int t = 0;

	if (file_export) {
		file = create_file();
		export_step(file, t);
	}

	const int id_x = id % q;
	const int id_y = id / q;

	double* line = (double*) malloc(6*(size_y/q)*sizeof(double));

	MPI_Request r[8] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
	MPI_Request s[8] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};

	for (t = 1; t < nb_steps; t++) {
		if (t == 1) {
			svdt = dt;
			dt = 0;
		}
		if (t == 2){
			dt = svdt / 2.;
		}

		// Travail par blocks, 1 block par processus

		// (1) Calcul bloc intÈrieur
		for (int j = id_y*(size_y/q)+1; j < (id_y+1)*(size_y/q)-1; j++)
		for (int i = id_x*(size_x/q)+1; i < (id_x+1)*(size_x/q)-1; i++)
			FORWARD(t, i, j);

		// (2) Attente rÈception bords
		MPI_Waitall(8, r, MPI_STATUSES_IGNORE);

		if (t != 1) {
			if (id_x != 0)
				for (int j = 0; j < size_y/q; j++)
					UPHY(t-1, id_x*(size_x/q)-1, id_y*(size_y/q)+j) = LINE(line, 3)[j];

			if (id_x != q-1) {
				for (int j = 0; j < size_y/q; j++)
					HPHY(t-1, (id_x+1)*(size_x/q), id_y*(size_y/q)+j) = LINE(line, 4)[j];

				for (int j = 0; j < size_y/q; j++)
					VPHY(t-1, (id_x+1)*(size_x/q), id_y*(size_y/q)+j) = LINE(line, 5)[j];
			}
		}

		// (3) Attente envoi anciens bords
		MPI_Waitall(8, s, MPI_STATUSES_IGNORE);

		// (4) RÈceptions suivantes
		if (id_y != q-1)
			MPI_Irecv(&VPHY(t, id_x*(size_x/q), (id_y+1)*(size_y/q)), size_x/q, MPI_DOUBLE, id_from_xy(id_x, id_y+1), 0, MPI_COMM_WORLD, r+0);

		if (id_y != 0) {
			MPI_Irecv(&HPHY(t, id_x*(size_x/q), id_y*(size_y/q)-1), size_x/q, MPI_DOUBLE, id_from_xy(id_x, id_y-1), 0, MPI_COMM_WORLD, r+1);
			MPI_Irecv(&UPHY(t, id_x*(size_x/q), id_y*(size_y/q)-1), size_x/q, MPI_DOUBLE, id_from_xy(id_x, id_y-1), 0, MPI_COMM_WORLD, r+2);
		}

		if (id_x != 0)
			MPI_Irecv(LINE(line, 3), size_y/q, MPI_DOUBLE, id_from_xy(id_x-1, id_y), 0, MPI_COMM_WORLD, r+3);

		if (id_x != q-1) {
			MPI_Irecv(LINE(line, 4), size_y/q, MPI_DOUBLE, id_from_xy(id_x+1, id_y), 0, MPI_COMM_WORLD, r+4);
			MPI_Irecv(LINE(line, 5), size_y/q, MPI_DOUBLE, id_from_xy(id_x+1, id_y), 0, MPI_COMM_WORLD, r+5);
		}

		if (id_y != 0 && id_x != 0)
			MPI_Irecv(&UPHY(t, id_x*(size_x/q)-1, id_y*(size_y/q))-1, 1, MPI_DOUBLE, id_from_xy(id_x-1, id_y-1), 0, MPI_COMM_WORLD, r+6);

		if (id_y != q-1 && id_x != q-1)
			MPI_Irecv(&VPHY(t, (id_x+1)*(size_x/q), (id_y+1)*(size_y/q)), 1, MPI_DOUBLE, id_from_xy(id_x+1, id_y+1), 0, MPI_COMM_WORLD, r+7);

		// (5) Calcul de mes bords
		for (int i = id_x*(size_x/q)+1; i < (id_x+1)*(size_x/q)-1; i++)
			FORWARD(t, i, id_y*(size_y/q));

		for (int i = id_x*(size_x/q)+1; i < (id_x+1)*(size_x/q)-1; i++)
			FORWARD(t, i, (id_y+1)*(size_y/q)-1);

		for (int j = id_y*(size_y/q)+1; j < (id_y+1)*(size_y/q); j++)
			FORWARD(t, id_x*(size_x/q), j);

		for (int j = id_y*(size_y/q); j < (id_y+1)*(size_y/q)-1; j++)
			FORWARD(t, (id_x+1)*(size_x/q)-1, j);

		FORWARD(t, id_x*(size_x/q), id_y*(size_y/q));
		FORWARD(t, (id_x+1)*(size_x/q)-1, (id_y+1)*(size_y/q)-1);

		// (6) Envoi de mes bords calculÈs
		if (id_y != q-1) {
			MPI_Isend(&HPHY(t, id_x*(size_x/q), (id_y+1)*(size_y/q)-1), size_x/q, MPI_DOUBLE, id_from_xy(id_x, id_y+1), 0, MPI_COMM_WORLD, s+0);
			MPI_Isend(&UPHY(t, id_x*(size_x/q), (id_y+1)*(size_y/q)-1), size_x/q, MPI_DOUBLE, id_from_xy(id_x, id_y+1), 0, MPI_COMM_WORLD, s+1);
		}
		if (id_y != 0)
			MPI_Isend(&VPHY(t, id_x*(size_x/q), id_y*(size_y/q)), size_x/q, MPI_DOUBLE, id_from_xy(id_x, id_y-1), 0, MPI_COMM_WORLD, s+2);

		if (id_x != 0) {
			for (int j = 0; j < size_y/q; j++)
				LINE(line, 0)[j] = HPHY(t, id_x*(size_x/q), id_y*(size_y/q)+j);
			MPI_Isend(LINE(line, 0), size_y/q, MPI_DOUBLE, id_from_xy(id_x-1, id_y), 0, MPI_COMM_WORLD, s+3);

			for (int j = 0; j < size_y/q; j++)
				LINE(line, 1)[j] = VPHY(t, id_x*(size_x/q), id_y*(size_y/q)+j);
			MPI_Isend(LINE(line, 1), size_y/q, MPI_DOUBLE, id_from_xy(id_x-1, id_y), 0, MPI_COMM_WORLD, s+4);
		}

		if (id_x != q-1) {
			for (int j = 0; j < size_y/q; j++)
				LINE(line, 2)[j] = UPHY(t, (id_x+1)*(size_x/q)-1, id_y*(size_y/q)+j);
			MPI_Isend(LINE(line, 2), size_y/q, MPI_DOUBLE, id_from_xy(id_x+1, id_y), 0, MPI_COMM_WORLD, s+5);
		}

		if (id_y != 0 && id_x != 0)
			MPI_Isend(&VPHY(t, id_x*(size_x/q), id_y*(size_y/q)), 1, MPI_DOUBLE, id_from_xy(id_x-1, id_y-1), 0, MPI_COMM_WORLD, s+6);

		if (id_y != q-1 && id_x != q-1)
			MPI_Isend(&UPHY(t, (id_x+1)*(size_x/q)-1, (id_y+1)*(size_y/q))-1, 1, MPI_DOUBLE, id_from_xy(id_x+1, id_y+1), 0, MPI_COMM_WORLD, s+7);


		if (file_export) {
			export_step(file, t);
		}

		if (t == 2) {
			dt = svdt;
		}
	}

	free(line);

	if (file_export) {
		finalize_export(file);
	}

	MPI_Barrier(MPI_COMM_WORLD);
}
