#include "export.h"
#include "shalw.h"
#include <math.h>
#include <stdio.h>

#include "simd_auto.h"
#include "simd_avx2.h"

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

inline void FORWARD(int t, int i, int j)
{
	HPHY(t, i, j) = hPhy_forward(t, i, j);
	HFIL(t, i, j) = hFil_forward(t, i, j);

	UPHY(t, i, j) = uPhy_forward(t, i, j);
	UFIL(t, i, j) = uFil_forward(t, i, j);

	VPHY(t, i, j) = vPhy_forward(t, i, j);
	VFIL(t, i, j) = vFil_forward(t, i, j);
}

void forward(void)
{
	FILE*  file = NULL;
	double svdt = 0.;
	int	t	= 0;

	if (file_export)
	{
		file = create_file();
		export_step(file, t);
	}

	for (t = 1; t < nb_steps; t++)
	{
		if (t == 1)
		{
			svdt = dt;
			dt   = 0;
		}
		if (t == 2)
		{
			dt = svdt / 2.;
		}


		if (t <= 2) {
			for (int j = 0; j < size_y; j++)
			for (int i = 0; i < size_x; i++)
				FORWARD(t, i, j);
		}
		else {
			// première ligne
			for (int i = 0; i < size_x; i++)
				FORWARD(t, i, 0);

			// rectangle intérieur
			for (int j = 1; j < size_y-1; j++)
			{
				FORWARD(t, 0, j);
				for (int i = 1; i < size_x-1; i++) {
					FORWARD_simd_auto(t, i, j);
				}
				FORWARD(t, size_x-1, j);

//				FORWARD(t, size_x-1, j);
				//				int i;
				//				for (i = 4; i+4-1 < size_x-1; i+=4) {
				//					FORWARD_simd_avx2(t, i, j);
				//				}
			}

			// dernière ligne
			for (int i = 0; i < size_x; i++)
				FORWARD(t, i, size_y-1);
		}

		if (file_export)
		{
			export_step(file, t);
		}

		if (t == 2)
		{
			dt = svdt;
		}
	}

	if (file_export)
	{
		finalize_export(file);
	}
}