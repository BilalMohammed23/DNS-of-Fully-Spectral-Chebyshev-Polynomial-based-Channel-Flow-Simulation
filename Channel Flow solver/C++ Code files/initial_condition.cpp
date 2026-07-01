#include "initial_condition.h"

#include <cmath>

FlowField create_initial_condition(const ChannelGrid& grid,
    double eps,
    double A)
{
    const double pi = 3.14159265358979323846;

    FlowField field;

    field.Nx = grid.Nx;
    field.Ny = grid.Nx;

    field.u.assign(field.Ny, std::vector<double>(field.Nx, 0.0));
    field.v.assign(field.Ny, std::vector<double>(field.Nx, 0.0));

    double alpha = 2.0 * pi / grid.L;

    for (int iy = 0; iy < field.Ny; iy++)
    {
        double Y = grid.y[iy];

        for (int ix = 0; ix < field.Nx; ix++)
        {
            double X = grid.x[ix];

            field.u[iy][ix] =
                (1.0 - Y * Y)
                - 4.0 * eps * A * Y * (1.0 - Y * Y) * std::sin(alpha * X);

            field.v[iy][ix] =
                -eps * A * alpha * std::pow((1.0 - Y * Y), 2.0) * std::cos(alpha * X);
        }
    }

    return field;
}