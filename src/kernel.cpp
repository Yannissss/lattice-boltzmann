#include <kernel.hpp>

#include <immintrin.h>

int accelerate_flow(const t_param params, t_speed cells, int *obstacles)
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
            (cells.speeds[3][ii + jj * params.nx] - w1) > 0.f &&
            (cells.speeds[6][ii + jj * params.nx] - w2) > 0.f &&
            (cells.speeds[7][ii + jj * params.nx] - w2) > 0.f)
        {
            /* increase 'east-side' densities */
            cells.speeds[1][ii + jj * params.nx] += w1;
            cells.speeds[5][ii + jj * params.nx] += w2;
            cells.speeds[8][ii + jj * params.nx] += w2;
            /* decrease 'west-side' densities */
            cells.speeds[3][ii + jj * params.nx] -= w1;
            cells.speeds[6][ii + jj * params.nx] -= w2;
            cells.speeds[7][ii + jj * params.nx] -= w2;
        }
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
                float local_density, acc[9];

                acc[0] += tmp_cells.speeds[0][ii + jj * params.nx];
                acc[1] += tmp_cells.speeds[1][ii + jj * params.nx];
                acc[2] += tmp_cells.speeds[2][ii + jj * params.nx];
                acc[3] += tmp_cells.speeds[3][ii + jj * params.nx];
                acc[4] += tmp_cells.speeds[4][ii + jj * params.nx];
                acc[5] += tmp_cells.speeds[5][ii + jj * params.nx];
                acc[6] += tmp_cells.speeds[6][ii + jj * params.nx];
                acc[7] += tmp_cells.speeds[7][ii + jj * params.nx];
                acc[8] += tmp_cells.speeds[8][ii + jj * params.nx];

                acc[0] = acc[0] + acc[1] + acc[2] + acc[3];
                acc[4] = acc[4] + acc[5] + acc[6] + acc[7];
                acc[0] = acc[0] + acc[4] + acc[8];

                local_density = acc[0];
                float inv_local_density = 1.f / local_density;

                /* compute x velocity component */
                float u_x = (tmp_cells.speeds[1][ii + jj * params.nx] +
                             tmp_cells.speeds[5][ii + jj * params.nx] +
                             tmp_cells.speeds[8][ii + jj * params.nx] -
                             (tmp_cells.speeds[3][ii + jj * params.nx] +
                              tmp_cells.speeds[6][ii + jj * params.nx] +
                              tmp_cells.speeds[7][ii + jj * params.nx])) *
                            inv_local_density;
                /* compute y velocity component */
                float u_y = (tmp_cells.speeds[2][ii + jj * params.nx] +
                             tmp_cells.speeds[5][ii + jj * params.nx] +
                             tmp_cells.speeds[6][ii + jj * params.nx] -
                             (tmp_cells.speeds[4][ii + jj * params.nx] +
                              tmp_cells.speeds[7][ii + jj * params.nx] +
                              tmp_cells.speeds[8][ii + jj * params.nx])) *
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

                /* relaxation step: unrolled now */
                cells.speeds[0][ii + jj * params.nx] =
                    tmp_cells.speeds[0][ii + jj * params.nx] +
                    params.omega *
                        (d_equ[0] - tmp_cells.speeds[0][ii + jj * params.nx]);

                cells.speeds[1][ii + jj * params.nx] =
                    tmp_cells.speeds[1][ii + jj * params.nx] +
                    params.omega *
                        (d_equ[1] - tmp_cells.speeds[1][ii + jj * params.nx]);

                cells.speeds[2][ii + jj * params.nx] =
                    tmp_cells.speeds[2][ii + jj * params.nx] +
                    params.omega *
                        (d_equ[2] - tmp_cells.speeds[2][ii + jj * params.nx]);

                cells.speeds[3][ii + jj * params.nx] =
                    tmp_cells.speeds[3][ii + jj * params.nx] +
                    params.omega *
                        (d_equ[3] - tmp_cells.speeds[3][ii + jj * params.nx]);

                cells.speeds[4][ii + jj * params.nx] =
                    tmp_cells.speeds[4][ii + jj * params.nx] +
                    params.omega *
                        (d_equ[4] - tmp_cells.speeds[4][ii + jj * params.nx]);

                cells.speeds[5][ii + jj * params.nx] =
                    tmp_cells.speeds[5][ii + jj * params.nx] +
                    params.omega *
                        (d_equ[5] - tmp_cells.speeds[5][ii + jj * params.nx]);

                cells.speeds[6][ii + jj * params.nx] =
                    tmp_cells.speeds[6][ii + jj * params.nx] +
                    params.omega *
                        (d_equ[6] - tmp_cells.speeds[6][ii + jj * params.nx]);

                cells.speeds[7][ii + jj * params.nx] =
                    tmp_cells.speeds[7][ii + jj * params.nx] +
                    params.omega *
                        (d_equ[7] - tmp_cells.speeds[7][ii + jj * params.nx]);

                cells.speeds[8][ii + jj * params.nx] =
                    tmp_cells.speeds[8][ii + jj * params.nx] +
                    params.omega *
                        (d_equ[8] - tmp_cells.speeds[8][ii + jj * params.nx]);
            }
        }
    }

    return EXIT_SUCCESS;
}

/* compute average velocity */
float av_velocity(const t_param params, t_speed cells, int *obstacles)
{
    /* initialise */
    __m512 ones = _mm512_set1_ps(1.0f);

    __m512 tot_u = _mm512_setzero_ps();
    __m512 tot_cells = _mm512_setzero_ps();

    /* loop over all non-blocked cells */
    for (int jj = 0; jj < params.ny; jj++)
    {
        __m512 local_speeds[NSPEEDS], u_x, u_y, u_u;
        __m512 local_density, inv_local_density;

        __m512 obsts_coeff;
        __m512i obsts;

        // Main computation
        for (int ii = 0; ii < params.nx; ii += 16)
        {
            /* ignore occupied cells */
            obsts = _mm512_load_si512(&obstacles[ii + jj * params.nx]);
            obsts_coeff = _mm512_cvtepi32_ps(obsts);

            /* local density total */

            local_speeds[0] =
                _mm512_load_ps(&cells.speeds[0][ii + jj * params.nx]);
            local_speeds[1] =
                _mm512_load_ps(&cells.speeds[1][ii + jj * params.nx]);
            local_speeds[2] =
                _mm512_load_ps(&cells.speeds[2][ii + jj * params.nx]);
            local_speeds[3] =
                _mm512_load_ps(&cells.speeds[3][ii + jj * params.nx]);
            local_speeds[4] =
                _mm512_load_ps(&cells.speeds[4][ii + jj * params.nx]);
            local_speeds[5] =
                _mm512_load_ps(&cells.speeds[5][ii + jj * params.nx]);
            local_speeds[6] =
                _mm512_load_ps(&cells.speeds[6][ii + jj * params.nx]);
            local_speeds[7] =
                _mm512_load_ps(&cells.speeds[7][ii + jj * params.nx]);
            local_speeds[8] =
                _mm512_load_ps(&cells.speeds[8][ii + jj * params.nx]);

            local_density = _mm512_setzero_ps();
            local_density = _mm512_add_ps(local_density, local_speeds[0]);
            local_density = _mm512_add_ps(local_density, local_speeds[1]);
            local_density = _mm512_add_ps(local_density, local_speeds[2]);
            local_density = _mm512_add_ps(local_density, local_speeds[3]);
            local_density = _mm512_add_ps(local_density, local_speeds[4]);
            local_density = _mm512_add_ps(local_density, local_speeds[5]);
            local_density = _mm512_add_ps(local_density, local_speeds[6]);
            local_density = _mm512_add_ps(local_density, local_speeds[7]);
            local_density = _mm512_add_ps(local_density, local_speeds[8]);

            inv_local_density = _mm512_div_ps(ones, local_density);

            /* x-component of velocity */
            u_x = _mm512_add_ps(u_x, local_speeds[1]);
            u_x = _mm512_add_ps(u_x, local_speeds[5]);
            u_x = _mm512_add_ps(u_x, local_speeds[8]);
            u_x = _mm512_sub_ps(u_x, local_speeds[3]);
            u_x = _mm512_sub_ps(u_x, local_speeds[6]);
            u_x = _mm512_sub_ps(u_x, local_speeds[7]);
            u_x = _mm512_mul_ps(u_x, inv_local_density);

            /* compute y velocity component */
            u_y = _mm512_add_ps(u_y, local_speeds[2]);
            u_y = _mm512_add_ps(u_y, local_speeds[5]);
            u_y = _mm512_add_ps(u_y, local_speeds[6]);
            u_y = _mm512_sub_ps(u_y, local_speeds[4]);
            u_y = _mm512_sub_ps(u_y, local_speeds[7]);
            u_y = _mm512_sub_ps(u_y, local_speeds[8]);
            u_y = _mm512_mul_ps(u_y, inv_local_density);

            // float u_u = sqrtf((u_x * u_x) + (u_y * u_y));
            u_u = _mm512_mul_ps(u_x, u_x);
            u_u = _mm512_add_ps(u_u, _mm512_mul_ps(u_y, u_y));
            u_u = _mm512_sqrt_ps(u_u);

            /* accumulate the norm of x- and y- velocity components */
            // tot_u = tot_u + obst_coeff * sqrtf((u_x * u_x) + (u_y * u_y));
            tot_u = _mm512_add_ps(tot_u, _mm512_mul_ps(obsts_coeff, u_u));

            /* increase counter of inspected cells */
            tot_cells = _mm512_add_ps(tot_cells, obsts_coeff);
        }
    }

    float u_acc, cells_acc;

    return _mm512_reduce_add_ps(tot_u) / _mm512_reduce_add_ps(tot_cells);
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
#pragma unroll
            for (int kk = 0; kk < NSPEEDS; kk++)
            {
                total += cells.speeds[kk][ii + jj * params.nx];
            }
        }
    }

    return total;
}
