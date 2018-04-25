#ifndef SIMD_AVX2_H
#define SIMD_AVX2_H

#include <immintrin.h>

#define ADD(x256, y256) (_mm256_add_pd((x256),(y256)))
#define SUB(x256, y256) (_mm256_sub_pd((x256),(y256)))
#define MUL(x256, y256) (_mm256_mul_pd((x256),(y256)))
#define DIV(x256, y256) (_mm256_div_pd((x256),(y256)))

// opère sur 4 doubles à la fois
inline void FORWARD_simd_avx2(int t, int i, int j)
{
	// m = moins, 0/1
	// p = plus, 0/1
	const __m256d UPHY_m1m1m1 = _mm256_loadu_pd(&UPHY(t - 1, i-1, j-1));

	const __m256d HPHY_m1p0p0 = _mm256_load_pd(&HPHY(t - 1, i, j));
	const __m256d UPHY_m1p0p0 = _mm256_load_pd(&UPHY(t - 1, i, j));
	const __m256d VPHY_m1p0p0 = _mm256_load_pd(&VPHY(t - 1, i, j));

	const __m256d UPHY_m1m1p0 = _mm256_loadu_pd(&UPHY(t - 1, i-1, j));

	const __m256d HFIL_m1p0p0 = _mm256_load_pd(&HFIL(t - 1, i, j));
	const __m256d UFIL_m1p0p0 = _mm256_load_pd(&UFIL(t - 1, i, j));
	const __m256d VFIL_m1p0p0 = _mm256_load_pd(&VFIL(t - 1, i, j));

	const __m256d VPHY_m1p0p1 = _mm256_loadu_pd(&VPHY(t - 1, i, j+1));

	const __m256d HPHY_m1p0m1 = _mm256_loadu_pd(&HPHY(t - 1, i, j-1));
	const __m256d UPHY_m1p0m1 = _mm256_loadu_pd(&UPHY(t - 1, i, j-1));

	const __m256d HPHY_m1p1p0 = _mm256_loadu_pd(&HPHY(t - 1, i+1, j));
	const __m256d VPHY_m1p1p0 = _mm256_loadu_pd(&VPHY(t - 1, i+1, j));

	const __m256d VPHY_m1p1p1 = _mm256_loadu_pd(&VPHY(t - 1, i+1, j+1));

	const __m256d HPHY_p0p0p0 = _mm256_load_pd(&HPHY(t, i, j));
	const __m256d UPHY_p0p0p0 = _mm256_load_pd(&UPHY(t, i, j));
	const __m256d VPHY_p0p0p0 = _mm256_load_pd(&VPHY(t, i, j));

	const __m256d _dt_hmoy = _mm256_set1_pd(dt*hmoy);
	const __m256d _dx = _mm256_set1_pd(dx);
	const __m256d _dy = _mm256_set1_pd(dy);

	// (UPHY(t - 1, i, j) - UPHY(t - 1, i - 1, j)) / dx
	const __m256d hy1 = DIV(SUB(UPHY_m1p0p0, UPHY_m1m1p0), _dx);

	// (VPHY(t - 1, i, j + 1) - VPHY(t - 1, i, j)) / dy
	const __m256d hy2 = DIV(SUB(VPHY_m1p0p1, VPHY_m1p0p0), _dy);

	//	HFIL(t - 1, i, j) -
	//				dt * hmoy * ((UPHY(t - 1, i, j) - UPHY(t - 1, i - 1, j)) / dx +
	//							 (VPHY(t - 1, i, j + 1) - VPHY(t - 1, i, j)) / dy)
	const __m256d hy = SUB(HFIL_m1p0p0, MUL(_dt_hmoy, ADD(hy1, hy2)));


	const __m256d _dt = _mm256_set1_pd(dt);
	const __m256d _grav_dx = _mm256_set1_pd(-grav / dx);
	const __m256d _pcor_4 = _mm256_set1_pd(pcor / 4.);
	const __m256d _dissip = _mm256_set1_pd(dissip);

	// (-grav / dx) * (HPHY(t - 1, i + 1, j) - HPHY(t - 1, i, j))
	const __m256d uy1 = MUL(_grav_dx, SUB(HPHY_m1p1p0, HPHY_m1p0p0));

	// (pcor / 4.) * (VPHY(t - 1, i, j) + VPHY(t - 1, i, j + 1) + VPHY(t - 1, i + 1, j) + VPHY(t - 1, i + 1, j + 1))
	const __m256d uy2 = MUL(_pcor_4, ADD(ADD(VPHY_m1p0p0, VPHY_m1p0p1), ADD(VPHY_m1p1p0, VPHY_m1p1p1)));

	// (dissip * UFIL(t - 1, i, j))
	const __m256d uy3 = MUL(_dissip, UFIL_m1p0p0);

	//	UFIL(t - 1, i, j) +
	//			   dt * ((-grav / dx) * (HPHY(t - 1, i + 1, j) - HPHY(t - 1, i, j)) +
	//					 (pcor / 4.) * (VPHY(t - 1, i, j) + VPHY(t - 1, i, j + 1) + VPHY(t - 1, i + 1, j) + VPHY(t - 1, i + 1, j + 1)) -
	//					 (dissip * UFIL(t - 1, i, j)))
	const __m256d uy = ADD(UFIL_m1p0p0, MUL(_dt, ADD(uy1, SUB(uy2, uy3))));


	// (-grav / dy) * (HPHY(t - 1, i, j) - HPHY(t - 1, i, j - 1))
	const __m256d vy1 = MUL(_grav_dx, SUB(HPHY_m1p0p0, HPHY_m1p0m1));

	// (pcor / 4.) * (UPHY(t - 1, i - 1, j - 1) + UPHY(t - 1, i - 1, j) + UPHY(t - 1, i, j - 1) + UPHY(t - 1, i, j))
	const __m256d vy2 = MUL(_pcor_4, ADD(ADD(UPHY_m1m1m1, UPHY_m1m1p0), ADD(UPHY_m1p0m1, UPHY_m1p0p0)));

	// (dissip * VFIL(t - 1, i, j))
	const __m256d vy3 = MUL(_dissip, VFIL_m1p0p0);

	//	VFIL(t - 1, i, j) +
	//			   dt * ((-grav / dy) * (HPHY(t - 1, i, j) - HPHY(t - 1, i, j - 1)) -
	//					 (pcor / 4.) * (UPHY(t - 1, i - 1, j - 1) + UPHY(t - 1, i - 1, j) + UPHY(t - 1, i, j - 1) + UPHY(t - 1, i, j)) -
	//					 (dissip * VFIL(t - 1, i, j)))
	const __m256d vy = ADD(VFIL_m1p0p0, MUL(_dt, SUB(vy1, SUB(vy2, vy3))));


	// HPHY(t, i, j) = hy;
	_mm256_store_pd(&HPHY(t, i, j), hy);

	// UPHY(t, i, j) = uy;
	_mm256_store_pd(&UPHY(t, i, j), uy);

	// VPHY(t, i, j) = vy;
	_mm256_store_pd(&VPHY(t, i, j), vy);

	const __m256d _a = _mm256_set1_pd(alpha);
	const __m256d _2 = _mm256_set1_pd(2);

	// HPHY(t - 1, i, j) + alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j))
	const __m256d hi = ADD(HPHY_m1p0p0, MUL(_a, SUB(HFIL_m1p0p0, ADD(MUL(_2, HPHY_m1p0p0), HPHY_p0p0p0))));

	// UPHY(t - 1, i, j) + alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j))
	const __m256d ui = ADD(UPHY_m1p0p0, MUL(_a, SUB(UFIL_m1p0p0, ADD(MUL(_2, UPHY_m1p0p0), UPHY_p0p0p0))));

	// VPHY(t - 1, i, j) + alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j))
	const __m256d vi = ADD(VPHY_m1p0p0, MUL(_a, SUB(VFIL_m1p0p0, ADD(MUL(_2, VPHY_m1p0p0), VPHY_p0p0p0))));

	// HFIL(t, i, j) = hi;
	_mm256_store_pd(&HFIL(t, i, j), hi);

	// UFIL(t, i, j) = ui;
	_mm256_store_pd(&UFIL(t, i, j), ui);

	// VFIL(t, i, j) = vi;
	_mm256_store_pd(&VFIL(t, i, j), vi);
}

#endif // SIMD_AVX2_H
