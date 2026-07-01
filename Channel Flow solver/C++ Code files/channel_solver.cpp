#include "channel_solver.h"

#include "timestep.h"
#include "nonlinear.h"
#include "viscous.h"
#include "pressure.h"
#include "rhs.h"
#include "uvp_solver.h"
#include "chebyshev.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

static double max_abs_real_matrix_solver(const RealMatrix& A)
{
    double max_value = 0.0;

    for (size_t iy = 0; iy < A.size(); iy++)
    {
        for (size_t ix = 0; ix < A[iy].size(); ix++)
        {
            max_value = std::max(max_value, std::abs(A[iy][ix]));
        }
    }

    return max_value;
}

static RealMatrix real_part_matrix_solver(const ComplexMatrix& input)
{
    int Ny = static_cast<int>(input.size());
    int Nx = static_cast<int>(input[0].size());

    RealMatrix output(Ny, std::vector<double>(Nx, 0.0));

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int ix = 0; ix < Nx; ix++)
        {
            output[iy][ix] = input[iy][ix].real();
        }
    }

    return output;
}

static ComplexMatrix compute_dudx_hatk_solver(const ChannelGrid& grid,
    const ComplexMatrix& u_hatk)
{
    int Ny = grid.Nx;
    int Nx = grid.Nx;

    ComplexMatrix dudx_hatk(Ny, std::vector<Complex>(Nx, Complex(0.0, 0.0)));

    Complex I(0.0, 1.0);

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int k = 0; k < Nx; k++)
        {
            dudx_hatk[iy][k] = I * grid.k_x[k] * u_hatk[iy][k];
        }
    }

    return dudx_hatk;
}

static ComplexMatrix compute_dvdy_hatk_solver(const ChannelGrid& grid,
    const ComplexMatrix& v_hatk)
{
    int Ny = grid.Nx;
    int Nx = grid.Nx;

    ComplexMatrix dvdy_hatk(Ny, std::vector<Complex>(Nx, Complex(0.0, 0.0)));

    for (int k = 0; k < Nx; k++)
    {
        std::vector<Complex> column(Ny);

        for (int iy = 0; iy < Ny; iy++)
        {
            column[iy] = v_hatk[iy][k];
        }

        std::vector<Complex> a = a_cheb_coeff(column, grid.theta);
        std::vector<Complex> b = b_cheb_coeff(a, Ny - 1);
        std::vector<Complex> derivative_column = cheb_eval_series(b, grid.y);

        for (int iy = 0; iy < Ny; iy++)
        {
            dvdy_hatk[iy][k] = derivative_column[iy];
        }
    }

    return dvdy_hatk;
}

static double compute_max_divergence_solver(const ChannelGrid& grid,
    const ComplexMatrix& u_hatk,
    const ComplexMatrix& v_hatk)
{
    ComplexMatrix dudx_hatk = compute_dudx_hatk_solver(grid, u_hatk);
    ComplexMatrix dvdy_hatk = compute_dvdy_hatk_solver(grid, v_hatk);

    RealMatrix dudx = real_part_matrix_solver(
        ifft_rows_complex_to_complex_fftw(dudx_hatk)
    );

    RealMatrix dvdy = real_part_matrix_solver(
        ifft_rows_complex_to_complex_fftw(dvdy_hatk)
    );

    int Ny = grid.Nx;
    int Nx = grid.Nx;

    double max_div = 0.0;

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int ix = 0; ix < Nx; ix++)
        {
            double div_value = dudx[iy][ix] + dvdy[iy][ix];
            max_div = std::max(max_div, std::abs(div_value));
        }
    }

    return max_div;
}

static FlowField make_flow_field_from_uv(const RealMatrix& u,
    const RealMatrix& v)
{
    FlowField field;

    field.Ny = static_cast<int>(u.size());
    field.Nx = static_cast<int>(u[0].size());
    field.u = u;
    field.v = v;

    return field;
}

ChannelSolverResult run_channel_solver(double final_time_target)
{
    const double pi = 3.14159265358979323846;

    int Nx = 64;
    double L = 2.0 * pi;

    ChannelGrid grid = create_channel_grid(Nx, L);

    double Re = 100;
    double K = -2.0 / Re;
    double gamma = 1.0 / Re;

    double eps = 0.01;
    double A = 1.0;

    FlowField field = create_initial_condition(grid, eps, A);

    ComplexMatrix u_hatk = fft_rows_real_to_complex_fftw(field.u);
    ComplexMatrix v_hatk = fft_rows_real_to_complex_fftw(field.v);

    double max_u = max_abs_real_matrix_solver(field.u);
    double max_v = max_abs_real_matrix_solver(field.v);

    NonlinearResult nonlinear_new = compute_nonlinear_terms(grid, u_hatk, v_hatk);
    ViscousResult viscous = compute_viscous_terms(grid, u_hatk, v_hatk);

    double max_div_0 = compute_max_divergence_solver(grid, u_hatk, v_hatk);

    PressureResult pressure0 = compute_pressure_0th_timestep(
        grid,
        nonlinear_new.H1,
        nonlinear_new.H2
    );

    NonlinearResult nonlinear_old = nonlinear_new;

    RealMatrix dPdx = pressure0.dPdx;
    RealMatrix dPdy = pressure0.dPdy;

    double t = 0.0;
    int count = 0;
    double max_divergence = max_div_0;

    std::cout << std::setprecision(12);
    std::cout << "Starting full channel solver" << std::endl;
    std::cout << "Nx = " << Nx << ", Re = " << Re << ", T = " << final_time_target << std::endl;
    std::cout << "Initial max divergence = " << max_div_0 << std::endl;
    std::cout << std::endl;

    while (t < final_time_target)
    {
        double dt = dt_channel(max_u, max_v, grid.dx, grid.dy_min, gamma);

        if (t + dt > final_time_target)
        {
            dt = final_time_target - t;
        }

        double beta = 2.0 / (gamma * dt);

        RHSResult rhs;

        if (count == 0)
        {
            rhs = compute_rhs_coefficients_initial(
                grid,
                K,
                gamma,
                dt,
                field,
                nonlinear_new,
                viscous,
                pressure0
            );
        }
        else
        {
            nonlinear_new = compute_nonlinear_terms(grid, u_hatk, v_hatk);
            viscous = compute_viscous_terms(grid, u_hatk, v_hatk);

            rhs = compute_rhs_coefficients_ab2(
                grid,
                K,
                gamma,
                dt,
                field,
                nonlinear_new,
                nonlinear_old,
                viscous,
                { RealMatrix(), dPdx, dPdy }
            );
        }

        UVPResult uvp = solve_uvp_chebyshev(
            grid,
            rhs.q_n_u,
            rhs.q_n_v,
            beta,
            gamma
        );

        u_hatk = uvp.u_hatk;
        v_hatk = uvp.v_hatk;

        field = make_flow_field_from_uv(uvp.u, uvp.v);

        dPdx = uvp.dPdx;
        dPdy = uvp.dPdy;

        max_u = max_abs_real_matrix_solver(field.u);
        max_v = max_abs_real_matrix_solver(field.v);

        max_divergence = compute_max_divergence_solver(grid, u_hatk, v_hatk);

        t += dt;
        count++;

        std::cout << std::setw(4) << count
            << ": t = " << std::setw(12) << t
            << ", dt = " << std::setw(12) << dt
            << ", max div = " << max_divergence
            << ", max |u| = " << max_u
            << ", max |v| = " << max_v
            << std::endl;

        nonlinear_old = nonlinear_new;
    }

    ChannelSolverResult result;

    result.grid = grid;

    result.u_final = field.u;
    result.v_final = field.v;
    result.u_hatk_final = u_hatk;
    result.v_hatk_final = v_hatk;
    result.dPdx_final = dPdx;
    result.dPdy_final = dPdy;
    result.final_time = t;
    result.step_count = count;
    result.max_divergence_final = max_divergence;

    return result;
}