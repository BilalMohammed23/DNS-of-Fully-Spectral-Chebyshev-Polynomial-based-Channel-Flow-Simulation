#pragma once

#include "grid.h"
#include "fourier.h"
#include "initial_condition.h"
#include "nonlinear.h"
#include "viscous.h"
#include "pressure.h"

struct RHSResult
{
    RealMatrix Q1;
    RealMatrix Q2;

    ComplexMatrix q_n_u;
    ComplexMatrix q_n_v;
};

RHSResult compute_rhs_coefficients_initial(const ChannelGrid& grid,
    double K,
    double gamma,
    double dt,
    const FlowField& field,
    const NonlinearResult& nonlinear,
    const ViscousResult& viscous,
    const PressureResult& pressure);

RHSResult compute_rhs_coefficients_ab2(const ChannelGrid& grid,
    double K,
    double gamma,
    double dt,
    const FlowField& field,
    const NonlinearResult& nonlinear_new,
    const NonlinearResult& nonlinear_old,
    const ViscousResult& viscous,
    const PressureResult& pressure);