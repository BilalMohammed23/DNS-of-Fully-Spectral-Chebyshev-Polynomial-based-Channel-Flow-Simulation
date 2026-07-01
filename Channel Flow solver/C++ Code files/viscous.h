#pragma once

#include "fourier.h"
#include "grid.h"

struct ViscousResult
{
    RealMatrix d2udx2;
    RealMatrix d2vdx2;
    RealMatrix d2udy2;
    RealMatrix d2vdy2;
};

ViscousResult compute_viscous_terms(const ChannelGrid& grid,
    const ComplexMatrix& uhatk,
    const ComplexMatrix& vhatk);