#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>
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

void forward(void) {
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

		// Calculs Blocs(t) -> Blocs(T+1)
		MPI_Request s[3] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};
		int computed = 0;
		while (computed != (U_COMPUTED|V_COMPUTED|H_COMPUTED)) {
			// Attendre d'avoir une réception
			int indx;
			MPI_Waitany(3, r, &indx, MPI_STATUS_IGNORE);

			//  mettre le flag recv_flags correspondant
			switch (indx) {
			case 0: received |= U_RECEIVED; r[0] = MPI_REQUEST_NULL; break;
			case 1: received |= V_RECEIVED; r[1] = MPI_REQUEST_NULL; break;
			case 2: received |= H_RECEIVED; r[2] = MPI_REQUEST_NULL; break;
			}

			// On regarde ce que l'on peut calculer
			if (received & (V_RECEIVED|H_RECEIVED) && (computed & U_COMPUTED) == 0) {
					// Travail par bandes, 1 bande par processus
					for (int j = id*(size_y/p); j < (id+1)*(size_y/p); j++) {
					for (int i = 0; i < size_x; i++) {
						UPHY(t, i, j) = uPhy_forward(t, i, j);}}

					// envoyer ma première ligne de UPHY
					if (id != p-1)
						MPI_Isend(&UPHY(t, 0, (id+1)*(size_y/p)-1), size_x, MPI_DOUBLE, id+1, U_TAG, MPI_COMM_WORLD, &s[0]);

					// recevoir leur dernière ligne de UPHY
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

					// envoyer ma première ligne de VPHY
					if (id != 0)
						MPI_Isend(&VPHY(t, 0, id*(size_y/p)), size_x, MPI_DOUBLE, id-1, V_TAG, MPI_COMM_WORLD, &s[1]);

					// recevoir leur dernière ligne de VPHY
					if (id != p-1)
						MPI_Irecv(&VPHY(t, 0, (id+1)*(size_y/p)), size_x, MPI_DOUBLE, id+1, V_TAG, MPI_COMM_WORLD, &r[1]);

					for (int j = id*(size_y/p); j < (id+1)*(size_y/p); j++) {
					for (int i = 0; i < size_x; i++) {
						VFIL(t, i, j) = vFil_forward(t, i, j);}}

					computed |= V_COMPUTED;
			}

			if (received & (U_RECEIVED|V_RECEIVED) && (computed & H_COMPUTED) == 0) {
					// Travail par bandes, 1 bande par processus
					for (int j = id*(size_y/p); j < (id+1)*(size_y/p); j++) {
					for (int i = 0; i < size_x; i++) {
						HPHY(t, i, j) = hPhy_forward(t, i, j);}}

					// envoyer ma première ligne de HPHY
					if (id != p-1)
						MPI_Isend(&HPHY(t, 0, (id+1)*(size_y/p)-1), size_x, MPI_DOUBLE, id+1, H_TAG, MPI_COMM_WORLD, &s[2]);

					// recevoir leur dernière ligne de HPHY
					if (id != 0)
						MPI_Irecv(&HPHY(t, 0, id*(size_y/p)-1), size_x, MPI_DOUBLE, id-1, H_TAG, MPI_COMM_WORLD, &r[2]);

					for (int j = id*(size_y/p); j < (id+1)*(size_y/p); j++) {
					for (int i = 0; i < size_x; i++) {
						HFIL(t, i, j) = hFil_forward(t, i, j);}}

					computed |= H_COMPUTED;
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
