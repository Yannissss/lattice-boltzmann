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
float av_velocity(const t_param params, t_speed *cells, int *obstacles);

/* Sum all the densities in the grid.
** The total should remain constant from one timestep to the next. */
float total_density(const t_param params, t_speed *cells);

/* calculate Reynolds number */
float calc_reynolds(const t_param params, t_speed *cells, int *obstacles);

#endif // KERNEL_HPP
