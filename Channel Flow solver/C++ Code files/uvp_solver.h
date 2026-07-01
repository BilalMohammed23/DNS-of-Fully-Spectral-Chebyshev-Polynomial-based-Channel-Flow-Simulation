#pragma once

#include "grid.h"
#include "fourier.h"

struct UVPResult
{
    ComplexMatrix u_hatk;
    ComplexMatrix v_hatk;
    ComplexMatrix P_hatk;

    RealMatrix u;
    RealMatrix v;

    RealMatrix dPdx;
    RealMatrix dPdy;
};

UVPResult solve_uvp_chebyshev(const ChannelGrid& grid,
    const ComplexMatrix& q_n_u,
    const ComplexMatrix& q_n_v,
    double beta,
    double gamma);