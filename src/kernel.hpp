#ifndef KERNEL_HPP
#define KERNEL_HPP

#include <cmath>
#include <cstdio>
#include <cstdlib>

// Speed resolution
#define NSPEEDS 9

/* struct to hold the parameter values */
typedef struct
{
    int nx;           /* no. of cells in x-direction */
    int ny;           /* no. of cells in y-direction */
    int maxIters;     /* no. of iterations */
    int reynolds_dim; /* dimension for Reynolds number */
    float density;    /* density per link */
    float accel;      /* density redistribution */
    float omega;      /* relaxation parameter */
} t_param;

/* struct to hold the 'speed' values */
typedef struct
{
    float speeds[NSPEEDS];
} t_speed;

/*
** The main calculation methods.
** timestep calls, in order, the functions:
** accelerate_flow(), propagate(), rebound() & collision()
*/
int accelerate_flow(const t_param params, t_speed *cells, int *obstacles);
int propagate(const t_param params, t_speed *cells, t_speed *tmp_cells);
int rebound(const t_param params, t_speed *cells, t_speed *tmp_cells,
            int *obstacles);
int collision(const t_param params, t_speed *cells, t_speed *tmp_cells,
              int *obstacles);
int write_values(const t_param params, t_speed *cells, int *obstacles,
                 float *av_vels);

static inline int timestep(const t_param params, t_speed *cells,
                           t_speed *tmp_cells, int *obstacles)
{
    accelerate_flow(params, cells, obstacles);
    propagate(params, cells, tmp_cells);
    rebound(params, cells, tmp_cells, obstacles);
    collision(params, cells, tmp_cells, obstacles);
    return EXIT_SUCCESS;
}

/* compute average velocity */
static inline float av_velocity(const t_param params, t_speed *cells,
                                int *obstacles)
{
    int tot_cells = 0; /* no. of cells used in calculation */
    float tot_u;       /* accumulated magnitudes of velocity for each cell */

    /* initialise */
    tot_u = 0.f;

    /* loop over all non-blocked cells */
    for (int jj = 0; jj < params.ny; jj++)
    {
        for (int ii = 0; ii < params.nx; ii++)
        {
            /* ignore occupied cells */
            if (!obstacles[ii + jj * params.nx])
            {
                /* local density total */
                float local_density = 0.f;

                local_density += cells[ii + jj * params.nx].speeds[0];
                local_density += cells[ii + jj * params.nx].speeds[1];
                local_density += cells[ii + jj * params.nx].speeds[2];
                local_density += cells[ii + jj * params.nx].speeds[3];
                local_density += cells[ii + jj * params.nx].speeds[4];
                local_density += cells[ii + jj * params.nx].speeds[5];
                local_density += cells[ii + jj * params.nx].speeds[6];
                local_density += cells[ii + jj * params.nx].speeds[7];
                local_density += cells[ii + jj * params.nx].speeds[8];

                float inv_local_density = 1.f / local_density;

                /* x-component of velocity */
                float u_x = (cells[ii + jj * params.nx].speeds[1] +
                             cells[ii + jj * params.nx].speeds[5] +
                             cells[ii + jj * params.nx].speeds[8] -
                             (cells[ii + jj * params.nx].speeds[3] +
                              cells[ii + jj * params.nx].speeds[6] +
                              cells[ii + jj * params.nx].speeds[7])) *
                            inv_local_density;
                /* compute y velocity component */
                float u_y = (cells[ii + jj * params.nx].speeds[2] +
                             cells[ii + jj * params.nx].speeds[5] +
                             cells[ii + jj * params.nx].speeds[6] -
                             (cells[ii + jj * params.nx].speeds[4] +
                              cells[ii + jj * params.nx].speeds[7] +
                              cells[ii + jj * params.nx].speeds[8])) *
                            inv_local_density;

                /* accumulate the norm of x- and y- velocity components */
                tot_u += sqrtf((u_x * u_x) + (u_y * u_y));
                /* increase counter of inspected cells */
                ++tot_cells;
            }
        }
    }

    return tot_u / (float)tot_cells;
}

/* Sum all the densities in the grid.
** The total should remain constant from one timestep to the next. */
float total_density(const t_param params, t_speed *cells);

/* calculate Reynolds number */
float calc_reynolds(const t_param params, t_speed *cells, int *obstacles);

#endif // KERNEL_HPP
