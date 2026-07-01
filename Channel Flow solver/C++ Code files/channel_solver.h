#pragma once

#include "grid.h"
#include "fourier.h"
#include "initial_condition.h"

struct ChannelSolverResult
{
    ChannelGrid grid;

    RealMatrix u_final;
    RealMatrix v_final;

    ComplexMatrix u_hatk_final;
    ComplexMatrix v_hatk_final;

    RealMatrix dPdx_final;
    RealMatrix dPdy_final;

    double final_time;
    int step_count;
    double max_divergence_final;
};

ChannelSolverResult run_channel_solver(double final_time_target);