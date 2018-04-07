#include <string>
extern double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
extern int	 size_x, size_y, size, nb_steps, size_block_x, size_block_y,
	size_block;
extern int		   start_block_x, start_block_y, end_block_x, end_block_y;
extern double	  dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
extern bool		   file_export;
extern std::string export_path;
extern int		   p, id, id_x, id_y, p_x, p_y;
extern bool		   async;
extern bool		   block;
extern int		   buffer_size;

// i:column, j:line
#define HFIL(t, i, j) hFil[(i) + (j)*size_x + ((t) % 2) * buffer_size]
#define UFIL(t, i, j) uFil[(i) + (j)*size_x + ((t) % 2) * buffer_size]
#define VFIL(t, i, j) vFil[(i) + (j)*size_x + ((t) % 2) * buffer_size]
#define HPHY(t, i, j) hPhy[(i) + (j)*size_x + ((t) % 2) * buffer_size]
#define UPHY(t, i, j) uPhy[(i) + (j)*size_x + ((t) % 2) * buffer_size]
#define VPHY(t, i, j) vPhy[(i) + (j)*size_x + ((t) % 2) * buffer_size]
