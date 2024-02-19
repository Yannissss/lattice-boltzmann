#include <kernel.hpp>

int timestep(const t_param params, t_speed cells, t_speed tmp_cells,
             int *obstacles)
{
    accelerate_flow(params, cells, obstacles);
    propagate(params, cells, tmp_cells);
    rebound(params, cells, tmp_cells, obstacles);
    collision(params, cells, tmp_cells, obstacles);
    return EXIT_SUCCESS;
}

int accelerate_flow(const t_param params, t_speed cells, int *obstacles)
{
    /* compute weighting factors */
    float w1 = params.density * params.accel / 9.f;
    float w2 = params.density * params.accel / 36.f;

    /* modify the 2nd row of the grid */
    const int jj = params.ny - 2;
    const int offset = jj * params.nx;

#pragma clang loop vectorize(enable)
    for (int ii = 0; ii < params.nx; ii++)
    {
        /* if the cell is not occupied and
        ** we don't send a negative density */
        int mask = (!obstacles[ii + offset]) &&
                   (cells.speeds[3][ii + offset] - w1) > 0.f &&
                   (cells.speeds[6][ii + offset] - w2) > 0.f &&
                   (cells.speeds[7][ii + offset] - w2) > 0.f;
        float maskf = float(mask);

        /* increase 'east-side' densities */
        cells.speeds[1][ii + offset] =
            w1 * maskf + cells.speeds[1][ii + offset];
        cells.speeds[5][ii + offset] =
            w2 * maskf + cells.speeds[5][ii + offset];
        cells.speeds[8][ii + offset] =
            w2 * maskf + cells.speeds[8][ii + offset];

        /* decrease 'west-side' densities */
        // cells.speeds[7][ii + offset] =
        //     w1 * float(mask) + cells.speeds[5][ii + offset];
        // cells.speeds[6][ii + offset] =
        //     w2 * float(-mask) + cells.speeds[6][ii + offset];
        // cells.speeds[7][ii + offset] =
        //     w2 * float(-mask) + cells.speeds[7][ii + offset];
    }

    return EXIT_SUCCESS;
}

int propagate(const t_param params, t_speed cells, t_speed tmp_cells)
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
            tmp_cells.speeds[0][ii + jj * params.nx] =
                cells.speeds[0][ii +
                                jj * params.nx]; /* central cell, no movement */
            tmp_cells.speeds[1][ii + jj * params.nx] =
                cells.speeds[1][x_w + jj * params.nx]; /* east */
            tmp_cells.speeds[2][ii + jj * params.nx] =
                cells.speeds[2][ii + y_s * params.nx]; /* north */
            tmp_cells.speeds[3][ii + jj * params.nx] =
                cells.speeds[3][x_e + jj * params.nx]; /* west */
            tmp_cells.speeds[4][ii + jj * params.nx] =
                cells.speeds[4][ii + y_n * params.nx]; /* south */
            tmp_cells.speeds[5][ii + jj * params.nx] =
                cells.speeds[5][x_w + y_s * params.nx]; /* north-east */
            tmp_cells.speeds[6][ii + jj * params.nx] =
                cells.speeds[6][x_e + y_s * params.nx]; /* north-west */
            tmp_cells.speeds[7][ii + jj * params.nx] =
                cells.speeds[7][x_e + y_n * params.nx]; /* south-west */
            tmp_cells.speeds[8][ii + jj * params.nx] =
                cells.speeds[8][x_w + y_n * params.nx]; /* south-east */
        }
    }

    return EXIT_SUCCESS;
}

int rebound(const t_param params, t_speed cells, t_speed tmp_cells,
            int *obstacles)
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
                cells.speeds[1][ii + jj * params.nx] =
                    tmp_cells.speeds[3][ii + jj * params.nx];
                cells.speeds[2][ii + jj * params.nx] =
                    tmp_cells.speeds[4][ii + jj * params.nx];
                cells.speeds[3][ii + jj * params.nx] =
                    tmp_cells.speeds[1][ii + jj * params.nx];
                cells.speeds[4][ii + jj * params.nx] =
                    tmp_cells.speeds[2][ii + jj * params.nx];
                cells.speeds[5][ii + jj * params.nx] =
                    tmp_cells.speeds[7][ii + jj * params.nx];
                cells.speeds[6][ii + jj * params.nx] =
                    tmp_cells.speeds[8][ii + jj * params.nx];
                cells.speeds[7][ii + jj * params.nx] =
                    tmp_cells.speeds[5][ii + jj * params.nx];
                cells.speeds[8][ii + jj * params.nx] =
                    tmp_cells.speeds[6][ii + jj * params.nx];
            }
        }
    }

    return EXIT_SUCCESS;
}

int collision(const t_param params, t_speed cells, t_speed tmp_cells,
              int *obstacles)
{
    const float c_sq = 1.f / 3.f; /* square of speed of sound */
    const float w0 = 4.f / 9.f;   /* weighting factor */
    const float w1 = 1.f / 9.f;   /* weighting factor */
    const float w2 = 1.f / 36.f;  /* weighting factor */

    /* loop over the cells in the grid
    ** NB the collision step is called after
    ** the propagate step and so values of interest
    ** are in the scratch-space grid */
    for (int jj = 0; jj < params.ny; jj++)
    {
        const int offset = jj * params.nx;

#pragma clang loop vectorize(enable)
        for (int ii = 0; ii < params.nx; ii++)
        {
            /* don't consider occupied cells */
            float obst_coeff = (1.f - obstacles[ii + jj * params.nx]);
            /* compute local density total */
            float local_density = 0.f;

            float speeds[NSPEEDS];

            // Unroll this loop since NSPEED is know at compile time
#pragma clang loop unroll(full)
            for (int kk = 0; kk < NSPEEDS; kk++)
            {
                speeds[kk] = tmp_cells.speeds[kk][ii + offset];
                local_density += speeds[kk];
            }
            float inv_local_density = 1.f / local_density;

            float u_x, u_y;

            /* x-component of velocity */
            u_x = speeds[1] + speeds[5] + speeds[8];
            u_x = inv_local_density * u_x - (speeds[3] + speeds[6] + speeds[7]);

            /* compute y velocity component */
            u_y = speeds[2] + speeds[5] + speeds[6];
            u_y = inv_local_density * u_x - (speeds[4] + speeds[7] + speeds[8]);

            /* velocity squared */
            float u_sq = u_x * u_x + (u_y * u_y);

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

            // Unroll this loop since NSPEED is know at compile time
#pragma clang loop unroll(full)
            for (int kk = 0; kk < NSPEEDS; kk++)
            {
                float aux;
                aux = d_equ[kk] - speeds[kk];
                aux = speeds[kk] + params.omega * aux;
                cells.speeds[kk][ii + offset] =
                    aux * obst_coeff +
                    (cells.speeds[kk][ii + offset] -
                     obst_coeff * cells.speeds[kk][ii + offset]);
            }
        }
    }

    return EXIT_SUCCESS;
}

float av_velocity(const t_param params, t_speed cells, int *obstacles)
{
    float tot_u;     /* accumulated magnitudes of velocity for each cell */
    float tot_cells; /* no. of cells used in calculation */

    /* initialise */
    tot_u = 0.f;
    tot_cells = 0.f;

    /* loop over all non-blocked cells */
    for (int jj = 0; jj < params.ny; jj++)
    {
        // Try to tell compiler to vectorize the ii loop
        const int offset = jj * params.nx;

#pragma clang loop vectorize(enable)
        for (int ii = 0; ii < params.nx; ii++)
        {
            float obst_coeff = (1.0 - obstacles[ii + offset]);

            /* ignore occupied cells */
            /* local density total */
            float local_density = 0.f;

            float speeds[NSPEEDS];

            // Unroll this loop since NSPEED is know at compile time
#pragma clang loop unroll(full)
            for (int kk = 0; kk < NSPEEDS; kk++)
            {
                speeds[kk] = cells.speeds[kk][ii + offset];
                local_density += speeds[kk];
            }
            float inv_local_density = 1.f / local_density;

            float u_x, u_y;

            /* x-component of velocity */
            u_x = speeds[1] + speeds[5] + speeds[8];
            u_x = inv_local_density * u_x - (speeds[3] + speeds[6] + speeds[7]);

            /* compute y velocity component */
            u_y = speeds[2] + speeds[5] + speeds[6];
            u_y = inv_local_density * u_x - (speeds[4] + speeds[7] + speeds[8]);

            /* accumulate the norm of x- and y- velocity components */
            u_x = u_x * u_x;
            u_y = u_y * u_y;

            tot_u = obst_coeff * u_x + tot_u;
            tot_u = obst_coeff * u_y + tot_u;

            /* increase counter of inspected cells */
            tot_cells = obst_coeff + tot_cells;
        }
    }

    return tot_u / tot_cells;
}

float calc_reynolds(const t_param params, t_speed cells, int *obstacles)
{
    const float viscosity = 1.f / 6.f * (2.f / params.omega - 1.f);

    return av_velocity(params, cells, obstacles) * params.reynolds_dim /
           viscosity;
}

float total_density(const t_param params, t_speed cells)
{
    float total = 0.f; /* accumulator */

    for (int jj = 0; jj < params.ny; jj++)
    {
        for (int ii = 0; ii < params.nx; ii++)
        {
#pragma clang loop unroll(full)
            for (int kk = 0; kk < NSPEEDS; kk++)
            {
                total += cells.speeds[kk][ii + jj * params.nx];
            }
        }
    }

    return total;
}
