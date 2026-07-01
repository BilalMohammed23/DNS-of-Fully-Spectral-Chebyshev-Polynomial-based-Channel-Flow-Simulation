#pragma once

#include "fourier.h"
#include "grid.h"

struct NonlinearResult
{
    RealMatrix H1;
    RealMatrix H2;

    ComplexMatrix dudx_hatk;
    ComplexMatrix dvdy_hatk;
};

NonlinearResult compute_nonlinear_terms(const ChannelGrid& grid,
    const ComplexMatrix& uhatk,
    const ComplexMatrix& vhatk);