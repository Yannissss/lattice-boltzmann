#include <kernel.hpp>

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
#pragma loop unroll
            for (int kk = 0; kk < NSPEEDS; kk++)
            {
                total += cells.speeds[kk][ii + jj * params.nx];
            }
        }
    }

    return total;
}
