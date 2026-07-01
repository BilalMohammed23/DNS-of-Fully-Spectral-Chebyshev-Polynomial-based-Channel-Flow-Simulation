#include "grid.h"

#include <cmath>

ChannelGrid create_channel_grid(int Nx, double L)
{
    const double pi = 3.14159265358979323846;

    ChannelGrid grid;

    grid.Nx = Nx;
    grid.L = L;
    grid.dx = L / static_cast<double>(Nx);

    grid.x.resize(Nx);
    grid.k_x.resize(Nx);
    grid.theta.resize(Nx);
    grid.y.resize(Nx);

    for (int j = 0; j < Nx; j++)
    {
        grid.x[j] = static_cast<double>(j) * grid.dx;
        grid.theta[j] = static_cast<double>(j) * pi / static_cast<double>(Nx - 1);
        grid.y[j] = std::cos(grid.theta[j]);
    }

    for (int j = 0; j < Nx; j++)
    {
        int l;

        if (j < Nx / 2)
        {
            l = j;
        }
        else
        {
            l = j - Nx;
        }

        grid.k_x[j] = (2.0 * pi / L) * static_cast<double>(l);
    }

    grid.dy_min = grid.y[0] - grid.y[1];

    return grid;
}