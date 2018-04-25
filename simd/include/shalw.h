#include <string>
extern double *__restrict__	hFil, *__restrict__ uFil, *__restrict__ vFil, *__restrict__ hPhy, *__restrict__ uPhy, *__restrict__ vPhy;
extern int		   size_x, size_y, nb_steps;
extern double	  dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
extern bool		   file_export;
extern std::string export_path;

#define HFIL(t, i, j) hFil[(i) + (j)*size_x + ((t) % 2) * size_x * size_y]
#define UFIL(t, i, j) uFil[(i) + (j)*size_x + ((t) % 2) * size_x * size_y]
#define VFIL(t, i, j) vFil[(i) + (j)*size_x + ((t) % 2) * size_x * size_y]
#define HPHY(t, i, j) hPhy[(i) + (j)*size_x + ((t) % 2) * size_x * size_y]
#define UPHY(t, i, j) uPhy[(i) + (j)*size_x + ((t) % 2) * size_x * size_y]
#define VPHY(t, i, j) vPhy[(i) + (j)*size_x + ((t) % 2) * size_x * size_y]
