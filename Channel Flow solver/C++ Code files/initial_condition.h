#pragma once

#include <vector>
#include "grid.h"

struct FlowField
{
    int Ny;
    int Nx;

    std::vector<std::vector<double>> u;
    std::vector<std::vector<double>> v;
};

FlowField create_initial_condition(const ChannelGrid& grid,
    double eps,
    double A);