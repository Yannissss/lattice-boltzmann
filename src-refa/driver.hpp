#ifndef DRIVER_HPP
#define DRIVER_HPP

#include <kernel.hpp>

void lattice_boltzmann(t_param const &params, t_speed *cells, t_speed *tmp_cells,
                       int *obstacles, float *av_vels);

#endif // DRIVER_HPP
