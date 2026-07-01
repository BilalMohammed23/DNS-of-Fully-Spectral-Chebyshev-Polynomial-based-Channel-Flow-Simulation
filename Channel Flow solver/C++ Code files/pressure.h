#pragma once

#include "fourier.h"
#include "grid.h"

struct PressureResult
{
    RealMatrix P;
    RealMatrix dPdx;
    RealMatrix dPdy;
};

PressureResult compute_pressure_0th_timestep(const ChannelGrid& grid,
    const RealMatrix& H1,
    const RealMatrix& H2);