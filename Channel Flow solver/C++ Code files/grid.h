#pragma once

#include <vector>

struct ChannelGrid
{
    int Nx;
    double L;
    double dx;
    double dy_min;

    std::vector<double> x;
    std::vector<double> k_x;
    std::vector<double> theta;
    std::vector<double> y;
};

ChannelGrid create_channel_grid(int Nx, double L);