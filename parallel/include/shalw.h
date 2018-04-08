#include <string>
#include <time.h>

extern double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
extern int	 size_x, size_y, size, nb_steps, size_block_x, size_block_y,
	size_block;
extern int		   start_block_x, start_block_y, end_block_x, end_block_y;
extern double	  dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
extern bool		   file_export;
extern int		   step_export;
extern std::string export_path;
extern int		   p, id, id_x, id_y, p_x, p_y;
extern bool		   async;
extern bool		   block;
extern int		   buffer_size;
extern clock_t	 start_time, end_time;

#define TIME(start, end) ((double)(end - start)) / CLOCKS_PER_SEC

#define ID0_TIME(message, instruction)                              \
	if (id == 0)                                                    \
	{                                                               \
		start_time = clock();                                       \
		(instruction);                                              \
		end_time = clock();                                         \
		printf(message " %.2fs\n",                                  \
			   ((double)(end_time - start_time)) / CLOCKS_PER_SEC); \
	}                                                               \
	else                                                            \
	{                                                               \
		(instruction);                                              \
	}

#define ID0(message, instruction) \
	if (id == 0)                  \
	{                             \
		printf(message);          \
		(instruction);            \
	}                             \
	else                          \
	{                             \
		(instruction);            \
	}

#define ID0_(instruction) \
	if (id == 0)          \
	{                     \
		(instruction);    \
	}

// i:column, j:line
#define HFIL(t, i, j) \
	hFil[(i) + (j) * (size_block_x + (block)*2) + ((t) % 2) * buffer_size]
#define UFIL(t, i, j) \
	uFil[(i) + (j) * (size_block_x + (block)*2) + ((t) % 2) * buffer_size]
#define VFIL(t, i, j) \
	vFil[(i) + (j) * (size_block_x + (block)*2) + ((t) % 2) * buffer_size]
#define HPHY(t, i, j) \
	hPhy[(i) + (j) * (size_block_x + (block)*2) + ((t) % 2) * buffer_size]
#define UPHY(t, i, j) \
	uPhy[(i) + (j) * (size_block_x + (block)*2) + ((t) % 2) * buffer_size]
#define VPHY(t, i, j) \
	vPhy[(i) + (j) * (size_block_x + (block)*2) + ((t) % 2) * buffer_size]
