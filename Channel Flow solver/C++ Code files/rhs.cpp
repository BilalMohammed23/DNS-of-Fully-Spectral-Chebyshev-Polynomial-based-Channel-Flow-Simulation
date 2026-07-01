#include "rhs.h"
#include "chebyshev.h"

#include <stdexcept>

static RealMatrix zero_real_matrix_rhs(int Ny, int Nx)
{
    return RealMatrix(Ny, std::vector<double>(Nx, 0.0));
}

static ComplexMatrix zero_complex_matrix_rhs(int Ny, int Nx)
{
    return ComplexMatrix(Ny, std::vector<Complex>(Nx, Complex(0.0, 0.0)));
}

static void check_real_matrix_size_rhs(const RealMatrix& A,
    int Ny,
    int Nx,
    const char* name)
{
    if (static_cast<int>(A.size()) != Ny)
    {
        throw std::runtime_error(std::string(name) + " has wrong number of rows.");
    }

    for (int iy = 0; iy < Ny; iy++)
    {
        if (static_cast<int>(A[iy].size()) != Nx)
        {
            throw std::runtime_error(std::string(name) + " has wrong number of columns.");
        }
    }
}

static ComplexMatrix chebyshev_coefficients_by_fourier_column(const ComplexMatrix& F_hatk,
    const std::vector<double>& theta)
{
    int Ny = static_cast<int>(F_hatk.size());
    int Nx = static_cast<int>(F_hatk[0].size());

    ComplexMatrix coeffs = zero_complex_matrix_rhs(Ny, Nx);

    for (int k = 0; k < Nx; k++)
    {
        std::vector<Complex> column(Ny);

        for (int iy = 0; iy < Ny; iy++)
        {
            column[iy] = F_hatk[iy][k];
        }

        std::vector<Complex> a = a_cheb_coeff(column, theta);

        for (int iy = 0; iy < Ny; iy++)
        {
            coeffs[iy][k] = a[iy];
        }
    }

    return coeffs;
}

RHSResult compute_rhs_coefficients_initial(const ChannelGrid& grid,
    double K,
    double gamma,
    double dt,
    const FlowField& field,
    const NonlinearResult& nonlinear,
    const ViscousResult& viscous,
    const PressureResult& pressure)
{
    int Nx = grid.Nx;
    int Ny = grid.Nx;

    check_real_matrix_size_rhs(field.u, Ny, Nx, "field.u");
    check_real_matrix_size_rhs(field.v, Ny, Nx, "field.v");

    RealMatrix Q1 = zero_real_matrix_rhs(Ny, Nx);
    RealMatrix Q2 = zero_real_matrix_rhs(Ny, Nx);

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int ix = 0; ix < Nx; ix++)
        {
            Q1[iy][ix] =
                -K
                - 0.5 * pressure.dPdx[iy][ix]
                + 0.5 * gamma * (viscous.d2udx2[iy][ix] + viscous.d2udy2[iy][ix])
                - nonlinear.H1[iy][ix]
                + field.u[iy][ix] / dt;

            Q2[iy][ix] =
                -0.5 * pressure.dPdy[iy][ix]
                + 0.5 * gamma * (viscous.d2vdx2[iy][ix] + viscous.d2vdy2[iy][ix])
                - nonlinear.H2[iy][ix]
                + field.v[iy][ix] / dt;
        }
    }

    ComplexMatrix Q1_hatk = fft_rows_real_to_complex_fftw(Q1);
    ComplexMatrix Q2_hatk = fft_rows_real_to_complex_fftw(Q2);

    RHSResult result;

    result.Q1 = Q1;
    result.Q2 = Q2;

    result.q_n_u = chebyshev_coefficients_by_fourier_column(Q1_hatk, grid.theta);
    result.q_n_v = chebyshev_coefficients_by_fourier_column(Q2_hatk, grid.theta);

    return result;
}

RHSResult compute_rhs_coefficients_ab2(const ChannelGrid& grid,
    double K,
    double gamma,
    double dt,
    const FlowField& field,
    const NonlinearResult& nonlinear_new,
    const NonlinearResult& nonlinear_old,
    const ViscousResult& viscous,
    const PressureResult& pressure)
{
    int Nx = grid.Nx;
    int Ny = grid.Nx;

    check_real_matrix_size_rhs(field.u, Ny, Nx, "field.u");
    check_real_matrix_size_rhs(field.v, Ny, Nx, "field.v");

    RealMatrix Q1 = zero_real_matrix_rhs(Ny, Nx);
    RealMatrix Q2 = zero_real_matrix_rhs(Ny, Nx);

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int ix = 0; ix < Nx; ix++)
        {
            Q1[iy][ix] =
                -K
                - 0.5 * pressure.dPdx[iy][ix]
                + 0.5 * gamma * (viscous.d2udx2[iy][ix] + viscous.d2udy2[iy][ix])
                - (1.5 * nonlinear_new.H1[iy][ix] - 0.5 * nonlinear_old.H1[iy][ix])
                + field.u[iy][ix] / dt;

            Q2[iy][ix] =
                -0.5 * pressure.dPdy[iy][ix]
                + 0.5 * gamma * (viscous.d2vdx2[iy][ix] + viscous.d2vdy2[iy][ix])
                - (1.5 * nonlinear_new.H2[iy][ix] - 0.5 * nonlinear_old.H2[iy][ix])
                + field.v[iy][ix] / dt;
        }
    }

    ComplexMatrix Q1_hatk = fft_rows_real_to_complex_fftw(Q1);
    ComplexMatrix Q2_hatk = fft_rows_real_to_complex_fftw(Q2);

    RHSResult result;

    result.Q1 = Q1;
    result.Q2 = Q2;

    result.q_n_u = chebyshev_coefficients_by_fourier_column(Q1_hatk, grid.theta);
    result.q_n_v = chebyshev_coefficients_by_fourier_column(Q2_hatk, grid.theta);

    return result;
}