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

static inline int accelerate_flow(const t_param params, t_speed *cells,
                                  int *obstacles)
{

    constexpr float ONE_9 = 1.f / 9.f;
    constexpr float ONE_36 = 1.f / 36.f;

    /* compute weighting factors */
    float w1 = params.density * params.accel * ONE_9;
    float w2 = params.density * params.accel * ONE_36;

    /* modify the 2nd row of the grid */
    int jj = params.ny - 2;

    for (int ii = 0; ii < params.nx; ii++)
    {
        /* if the cell is not occupied and
        ** we don't send a negative density */
        if (!obstacles[ii + jj * params.nx] &&
            (cells[ii + jj * params.nx].speeds[3] - w1) > 0.f &&
            (cells[ii + jj * params.nx].speeds[6] - w2) > 0.f &&
            (cells[ii + jj * params.nx].speeds[7] - w2) > 0.f)
        {
            /* increase 'east-side' densities */
            cells[ii + jj * params.nx].speeds[1] += w1;
            cells[ii + jj * params.nx].speeds[5] += w2;
            cells[ii + jj * params.nx].speeds[8] += w2;
            /* decrease 'west-side' densities */
            cells[ii + jj * params.nx].speeds[3] -= w1;
            cells[ii + jj * params.nx].speeds[6] -= w2;
            cells[ii + jj * params.nx].speeds[7] -= w2;
        }
    }

    return EXIT_SUCCESS;
}

static inline int propagate(const t_param params, t_speed *cells,
                            t_speed *tmp_cells)
{
    /* loop over _all_ cells */
    for (int jj = 0; jj < params.ny; jj++)
    {
        for (int ii = 0; ii < params.nx; ii++)
        {
            /* determine indices of axis-direction neighbours
            ** respecting periodic boundary conditions (wrap around) */
            int y_n = (jj + 1) % params.ny;
            int x_e = (ii + 1) % params.nx;
            int y_s = (jj == 0) ? (jj + params.ny - 1) : (jj - 1);
            int x_w = (ii == 0) ? (ii + params.nx - 1) : (ii - 1);
            /* propagate densities from neighbouring cells, following
            ** appropriate directions of travel and writing into
            ** scratch space grid */
            tmp_cells[ii + jj * params.nx].speeds[0] =
                cells[ii + jj * params.nx]
                    .speeds[0]; /* central cell, no movement */
            tmp_cells[ii + jj * params.nx].speeds[1] =
                cells[x_w + jj * params.nx].speeds[1]; /* east */
            tmp_cells[ii + jj * params.nx].speeds[2] =
                cells[ii + y_s * params.nx].speeds[2]; /* north */
            tmp_cells[ii + jj * params.nx].speeds[3] =
                cells[x_e + jj * params.nx].speeds[3]; /* west */
            tmp_cells[ii + jj * params.nx].speeds[4] =
                cells[ii + y_n * params.nx].speeds[4]; /* south */
            tmp_cells[ii + jj * params.nx].speeds[5] =
                cells[x_w + y_s * params.nx].speeds[5]; /* north-east */
            tmp_cells[ii + jj * params.nx].speeds[6] =
                cells[x_e + y_s * params.nx].speeds[6]; /* north-west */
            tmp_cells[ii + jj * params.nx].speeds[7] =
                cells[x_e + y_n * params.nx].speeds[7]; /* south-west */
            tmp_cells[ii + jj * params.nx].speeds[8] =
                cells[x_w + y_n * params.nx].speeds[8]; /* south-east */
        }
    }

    return EXIT_SUCCESS;
}

static inline int rebound(const t_param params, t_speed *cells,
                          t_speed *tmp_cells, int *obstacles)
{
    /* loop over the cells in the grid */
    for (int jj = 0; jj < params.ny; jj++)
    {
        for (int ii = 0; ii < params.nx; ii++)
        {
            /* if the cell contains an obstacle */
            if (obstacles[jj * params.nx + ii])
            {
                /* called after propagate, so taking values from scratch space
                ** mirroring, and writing into main grid */
                cells[ii + jj * params.nx].speeds[1] =
                    tmp_cells[ii + jj * params.nx].speeds[3];
                cells[ii + jj * params.nx].speeds[2] =
                    tmp_cells[ii + jj * params.nx].speeds[4];
                cells[ii + jj * params.nx].speeds[3] =
                    tmp_cells[ii + jj * params.nx].speeds[1];
                cells[ii + jj * params.nx].speeds[4] =
                    tmp_cells[ii + jj * params.nx].speeds[2];
                cells[ii + jj * params.nx].speeds[5] =
                    tmp_cells[ii + jj * params.nx].speeds[7];
                cells[ii + jj * params.nx].speeds[6] =
                    tmp_cells[ii + jj * params.nx].speeds[8];
                cells[ii + jj * params.nx].speeds[7] =
                    tmp_cells[ii + jj * params.nx].speeds[5];
                cells[ii + jj * params.nx].speeds[8] =
                    tmp_cells[ii + jj * params.nx].speeds[6];
            }
        }
    }

    return EXIT_SUCCESS;
}

static inline int collision(const t_param params, t_speed *cells,
                            t_speed *tmp_cells, int *obstacles)
{
    constexpr float c_sq = 1.f / 3.f; /* square of speed of sound */
    constexpr float w0 = 4.f / 9.f;   /* weighting factor */
    constexpr float w1 = 1.f / 9.f;   /* weighting factor */
    constexpr float w2 = 1.f / 36.f;  /* weighting factor */

    /* loop over the cells in the grid
    ** NB the collision step is called after
    ** the propagate step and so values of interest
    ** are in the scratch-space grid */
    for (int jj = 0; jj < params.ny; jj++)
    {
        for (int ii = 0; ii < params.nx; ii++)
        {
            /* don't consider occupied cells */
            if (!obstacles[ii + jj * params.nx])
            {
                /* compute local density total */
                float local_density = 0.f;

                local_density += tmp_cells[ii + jj * params.nx].speeds[0];
                local_density += tmp_cells[ii + jj * params.nx].speeds[1];
                local_density += tmp_cells[ii + jj * params.nx].speeds[2];
                local_density += tmp_cells[ii + jj * params.nx].speeds[3];
                local_density += tmp_cells[ii + jj * params.nx].speeds[4];
                local_density += tmp_cells[ii + jj * params.nx].speeds[5];
                local_density += tmp_cells[ii + jj * params.nx].speeds[6];
                local_density += tmp_cells[ii + jj * params.nx].speeds[7];
                local_density += tmp_cells[ii + jj * params.nx].speeds[8];

                float inv_local_density = 1.f / local_density;

                /* compute x velocity component */
                float u_x = (tmp_cells[ii + jj * params.nx].speeds[1] +
                             tmp_cells[ii + jj * params.nx].speeds[5] +
                             tmp_cells[ii + jj * params.nx].speeds[8] -
                             (tmp_cells[ii + jj * params.nx].speeds[3] +
                              tmp_cells[ii + jj * params.nx].speeds[6] +
                              tmp_cells[ii + jj * params.nx].speeds[7])) *
                            inv_local_density;
                /* compute y velocity component */
                float u_y = (tmp_cells[ii + jj * params.nx].speeds[2] +
                             tmp_cells[ii + jj * params.nx].speeds[5] +
                             tmp_cells[ii + jj * params.nx].speeds[6] -
                             (tmp_cells[ii + jj * params.nx].speeds[4] +
                              tmp_cells[ii + jj * params.nx].speeds[7] +
                              tmp_cells[ii + jj * params.nx].speeds[8])) *
                            inv_local_density;

                /* velocity squared */
                float u_sq = u_x * u_x + u_y * u_y;

                /* directional velocity components */
                float u[NSPEEDS];
                u[1] = u_x;        /* east */
                u[2] = u_y;        /* north */
                u[3] = -u_x;       /* west */
                u[4] = -u_y;       /* south */
                u[5] = u_x + u_y;  /* north-east */
                u[6] = -u_x + u_y; /* north-west */
                u[7] = -u_x - u_y; /* south-west */
                u[8] = u_x - u_y;  /* south-east */

                /* equilibrium densities */
                float d_equ[NSPEEDS];
                /* zero velocity density: weight w0 */
                d_equ[0] = w0 * local_density * (1.f - u_sq / (2.f * c_sq));
                /* axis speeds: weight w1 */
                d_equ[1] =
                    w1 * local_density *
                    (1.f + u[1] / c_sq + (u[1] * u[1]) / (2.f * c_sq * c_sq) -
                     u_sq / (2.f * c_sq));
                d_equ[2] =
                    w1 * local_density *
                    (1.f + u[2] / c_sq + (u[2] * u[2]) / (2.f * c_sq * c_sq) -
                     u_sq / (2.f * c_sq));
                d_equ[3] =
                    w1 * local_density *
                    (1.f + u[3] / c_sq + (u[3] * u[3]) / (2.f * c_sq * c_sq) -
                     u_sq / (2.f * c_sq));
                d_equ[4] =
                    w1 * local_density *
                    (1.f + u[4] / c_sq + (u[4] * u[4]) / (2.f * c_sq * c_sq) -
                     u_sq / (2.f * c_sq));
                /* diagonal speeds: weight w2 */
                d_equ[5] =
                    w2 * local_density *
                    (1.f + u[5] / c_sq + (u[5] * u[5]) / (2.f * c_sq * c_sq) -
                     u_sq / (2.f * c_sq));
                d_equ[6] =
                    w2 * local_density *
                    (1.f + u[6] / c_sq + (u[6] * u[6]) / (2.f * c_sq * c_sq) -
                     u_sq / (2.f * c_sq));
                d_equ[7] =
                    w2 * local_density *
                    (1.f + u[7] / c_sq + (u[7] * u[7]) / (2.f * c_sq * c_sq) -
                     u_sq / (2.f * c_sq));
                d_equ[8] =
                    w2 * local_density *
                    (1.f + u[8] / c_sq + (u[8] * u[8]) / (2.f * c_sq * c_sq) -
                     u_sq / (2.f * c_sq));

                /* relaxation step */
                for (int kk = 0; kk < NSPEEDS; kk++)
                {
                    cells[ii + jj * params.nx].speeds[kk] =
                        tmp_cells[ii + jj * params.nx].speeds[kk] +
                        params.omega *
                            (d_equ[kk] -
                             tmp_cells[ii + jj * params.nx].speeds[kk]);
                }
            }
        }
    }

    return EXIT_SUCCESS;
}

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
