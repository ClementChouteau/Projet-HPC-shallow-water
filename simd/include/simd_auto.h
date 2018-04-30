#ifndef SIMD_AUTO_H
#define SIMD_AUTO_H

// t > 2
inline double hFil_forward_simd_auto(int t, int i, int j)
{
	return HPHY(t - 1, i, j) + alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

// t > 2
inline double uFil_forward_simd_auto(int t, int i, int j)
{
	return UPHY(t - 1, i, j) + alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

// t > 2
inline double vFil_forward_simd_auto(int t, int i, int j)
{
	return VPHY(t - 1, i, j) + alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
}

// i > 0 et j < size_y - 1
inline double hPhy_forward_simd_auto(int t, int i, int j)
{
	const double c = UPHY(t - 1, i - 1, j);
	const double d = VPHY(t - 1, i, j + 1);

	return HFIL(t - 1, i, j) - dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx + (d - VPHY(t - 1, i, j)) / dy);
}

// i < size_x - 1 et j < size_y - 1
inline double uPhy_forward_simd_auto(int t, int i, int j)
{
	const double b = HPHY(t - 1, i + 1, j);
	const double e = VPHY(t - 1, i, j + 1);
	const double f = VPHY(t - 1, i + 1, j);
	const double g = VPHY(t - 1, i + 1, j + 1);

	return UFIL(t - 1, i, j) + dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) + (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) - (dissip * UFIL(t - 1, i, j)));
}

// j > 0 et i > 0
inline double vPhy_forward_simd_auto(int t, int i, int j)
{
	const double c = HPHY(t - 1, i, j - 1);
	const double d = UPHY(t - 1, i - 1, j - 1);
	const double e = UPHY(t - 1, i - 1, j);
	const double f = UPHY(t - 1, i, j - 1);

	return VFIL(t - 1, i, j) + dt * ((-grav / dy) * (HPHY(t - 1, i, j) - c) - (pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) - (dissip * VFIL(t - 1, i, j)));
}

inline void FORWARD_simd_auto(int t, int i, int j)
{
	const double uy = uPhy_forward_simd_auto(t, i, j);
	const double vy = vPhy_forward_simd_auto(t, i, j);
	const double hy = hPhy_forward_simd_auto(t, i, j);

	UPHY(t, i, j) = uy;
	VPHY(t, i, j) = vy;
	HPHY(t, i, j) = hy;

	const double ui = uFil_forward_simd_auto(t, i, j);
	const double vi = vFil_forward_simd_auto(t, i, j);
	const double hi = hFil_forward_simd_auto(t, i, j);

	UFIL(t, i, j) = ui;
	VFIL(t, i, j) = vi;
	HFIL(t, i, j) = hi;
}

#endif // SIMD_AUTO_H
