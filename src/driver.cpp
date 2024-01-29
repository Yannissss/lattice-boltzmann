#include <driver.hpp>

void lattice_boltzmann(t_param const &params, t_speed cells, t_speed tmp_cells,
                       int *obstacles, float *av_vels)
{
    for (int tt = 0; tt < params.maxIters; tt++)
    {
        timestep(params, cells, tmp_cells, obstacles);
        av_vels[tt] = av_velocity(params, cells, obstacles);
#ifdef DEBUG
        printf("==timestep: %d==\n", tt);
        printf("av velocity: %.12E\n", av_vels[tt]);
        printf("tot density: %.12E\n", total_density(params, cells));
#endif
        /** Echange des données */
    }
}