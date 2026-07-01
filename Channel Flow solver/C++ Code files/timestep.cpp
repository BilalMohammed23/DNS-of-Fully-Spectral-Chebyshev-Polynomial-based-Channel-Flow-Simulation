#include "timestep.h"

#include <algorithm>

double dt_channel(double max_u,
    double max_v,
    double dx,
    double dy_min,
    double gamma)
{
    const double eps0 = 1.0e-12;
    const double CFL = 0.2;

    double dt_x = CFL * dx / std::max(max_u, eps0);
    double dt_y = CFL * dy_min / std::max(max_v, eps0);
    double dt_conv = std::min(dt_x, dt_y);

    double dt_diff = 0.5 / (gamma * (1.0 / (dx * dx) + 1.0 / (dy_min * dy_min)));

    return std::min(dt_conv, dt_diff);
}