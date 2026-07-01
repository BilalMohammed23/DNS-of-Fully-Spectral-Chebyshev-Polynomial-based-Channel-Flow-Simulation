#include "uvp_solver.h"
#include "chebyshev.h"
#include "linear_solver.h"

#include <stdexcept>
#include <cmath>

static ComplexMatrix zero_complex_matrix_uvp(int Ny, int Nx)
{
    return ComplexMatrix(Ny, std::vector<Complex>(Nx, Complex(0.0, 0.0)));
}

static RealMatrix zero_real_matrix_uvp(int Ny, int Nx)
{
    return RealMatrix(Ny, std::vector<double>(Nx, 0.0));
}

static RealMatrix real_part_matrix_uvp(const ComplexMatrix& input)
{
    int Ny = static_cast<int>(input.size());
    int Nx = static_cast<int>(input[0].size());

    RealMatrix output = zero_real_matrix_uvp(Ny, Nx);

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int ix = 0; ix < Nx; ix++)
        {
            output[iy][ix] = input[iy][ix].real();
        }
    }

    return output;
}

static ComplexVector get_fourier_column(const ComplexMatrix& A,
    int k)
{
    int Ny = static_cast<int>(A.size());

    ComplexVector column(Ny, Complex(0.0, 0.0));

    for (int iy = 0; iy < Ny; iy++)
    {
        column[iy] = A[iy][k];
    }

    return column;
}

static Complex q_at(const ComplexVector& q,
    int matlab_index)
{
    return q[matlab_index - 1];
}

static ComplexVector solve_global_matrix_for_mode(const ChannelGrid& grid,
    int ix,
    const ComplexMatrix& q_n_u,
    const ComplexMatrix& q_n_v,
    double beta,
    double gamma,
    bool use_k0_matrix)
{
    int Nx = grid.Nx;
    int Ntot = 3 * Nx;

    ComplexDenseMatrix A(Ntot, ComplexVector(Ntot, Complex(0.0, 0.0)));
    ComplexVector b(Ntot, Complex(0.0, 0.0));

    double kx_current = grid.k_x[ix];
    double alpha_eqn = kx_current * kx_current + beta;

    ComplexVector q_n_mu = get_fourier_column(q_n_u, ix);
    ComplexVector q_n_mv = get_fourier_column(q_n_v, ix);

    Complex I(0.0, 1.0);

    double c_n = 1.0;
    double c_n_1 = 1.0;

    for (int i = 1; i <= Ntot; i++)
    {
        int i0 = i - 1;

        // u-RHS block
        if ((i > 2) && (i < Nx + 1))
        {
            double c_n_2 = 1.0;

            if (i - 3 == 0)
            {
                c_n_2 = 2.0;
            }

            Complex RHS1_u =
                -((c_n_1 * 2.0) / (2.0 * gamma * (i - 2)))
                * (c_n_2 * q_at(q_n_mu, i - 2) - q_at(q_n_mu, i));

            Complex RHS2_u;

            if ((i == Nx) || (i == Nx - 1))
            {
                RHS2_u =
                    (2.0 / (2.0 * gamma * i))
                    * (c_n * q_at(q_n_mu, i));
            }
            else
            {
                RHS2_u =
                    (2.0 / (2.0 * gamma * i))
                    * (c_n * q_at(q_n_mu, i) - q_at(q_n_mu, i + 2));
            }

            b[i0] = RHS1_u + RHS2_u;
        }

        // v-RHS block
        if ((i > Nx + 2) && (i < 2 * Nx + 1))
        {
            double c_n_2 = 1.0;

            if (i - (Nx + 1) - 2 == 0)
            {
                c_n_2 = 2.0;
            }

            Complex RHS1_v =
                (c_n_1 * q_at(q_n_mv, i - Nx)
                    - c_n_1 * c_n_2 * q_at(q_n_mv, i - Nx - 2))
                / ((i - Nx - 2) * gamma);

            Complex RHS2_v;

            if ((i == 2 * Nx) || (i == 2 * Nx - 1))
            {
                RHS2_v =
                    -(-c_n * q_at(q_n_mv, i - Nx))
                    / ((i - Nx) * gamma);
            }
            else
            {
                RHS2_v =
                    -(q_at(q_n_mv, i - Nx + 2) - c_n * q_at(q_n_mv, i - Nx))
                    / ((i - Nx) * gamma);
            }

            b[i0] = RHS1_v + RHS2_v;
        }

        for (int j = 1; j <= Ntot; j++)
        {
            int j0 = j - 1;

            // ------------------------------------------------------------
            // u-velocity block
            // ------------------------------------------------------------

            // u boundary rows
            if ((i == 1) && (j < Nx + 1))
            {
                A[i0][j0] = Complex(1.0, 0.0);
            }
            else if ((i == 2) && (j < Nx + 1))
            {
                int exponent = j - 1;
                double value = (exponent % 2 == 0) ? 1.0 : -1.0;
                A[i0][j0] = Complex(value, 0.0);
            }

            // u interior rows
            if ((i > 2) && (i < Nx + 1) && (j < Nx + 1))
            {
                double c_n_2 = 1.0;

                if (i - 3 == 0)
                {
                    c_n_2 = 2.0;
                }

                if (i == j)
                {
                    A[i0][j0] =
                        (c_n_1 * alpha_eqn / (2.0 * (i - 2)))
                        + (alpha_eqn * c_n / (2.0 * i))
                        + 2.0 * (i - 1);
                }

                if ((j < i) && (j == i - 2))
                {
                    A[i0][j0] =
                        -(alpha_eqn * c_n_1 * c_n_2 / (2.0 * (i - 2)));
                }
                else if ((j > i) && (j == i + 2))
                {
                    A[i0][j0] =
                        -alpha_eqn / (2.0 * i);
                }
            }

            // ------------------------------------------------------------
            // u-pressure block
            // This block is removed in the special k = 0 matrix.
            // ------------------------------------------------------------
            if (!use_k0_matrix)
            {
                if ((i > 2) && (i < Nx + 1) && (j > 2 * Nx))
                {
                    double c_n_2 = 1.0;

                    if (i - 3 == 0)
                    {
                        c_n_2 = 2.0;
                    }

                    if (i == (j - 2 * Nx))
                    {
                        A[i0][j0] =
                            ((c_n_1 * I * kx_current) / (2.0 * gamma * (i - 2)))
                            + ((c_n * I * kx_current) / (2.0 * gamma * i));
                    }

                    if ((j < (i + 2 * Nx)) && (j == (i + 2 * Nx) - 2))
                    {
                        A[i0][j0] =
                            -(c_n_1 * I * kx_current * c_n_2)
                            / (2.0 * gamma * (i - 2));
                    }
                    else if ((j > (i + 2 * Nx)) && (j == (i + 2 * Nx) + 2))
                    {
                        A[i0][j0] =
                            -(I * kx_current) / (2.0 * gamma * i);
                    }
                }
            }

            // ------------------------------------------------------------
            // v-velocity block
            // ------------------------------------------------------------

            // v boundary rows
            if ((i == Nx + 1) && (j > Nx) && (j < 2 * Nx + 1))
            {
                A[i0][j0] = Complex(1.0, 0.0);
            }
            else if ((i == Nx + 2) && (j > Nx) && (j < 2 * Nx + 1))
            {
                int exponent = j - Nx - 1;
                double value = (exponent % 2 == 0) ? 1.0 : -1.0;
                A[i0][j0] = Complex(value, 0.0);
            }

            // v interior rows
            if ((i > Nx + 2)
                && (i < 2 * Nx + 1)
                && (j > Nx)
                && (j < 2 * Nx + 1))
            {
                double c_n_2 = 1.0;

                if (i - (Nx + 1) - 2 == 0)
                {
                    c_n_2 = 2.0;
                }

                if (i == j)
                {
                    A[i0][j0] =
                        (c_n_1 * alpha_eqn / (2.0 * (i - Nx - 2)))
                        + (alpha_eqn * c_n / (2.0 * (i - Nx)))
                        + 2.0 * (i - Nx - 1);
                }

                if ((j < i) && (j == i - 2))
                {
                    A[i0][j0] =
                        -(alpha_eqn * c_n_1 * c_n_2 / (2.0 * (i - Nx - 2)));
                }
                else if ((j > i) && (j == i + 2))
                {
                    A[i0][j0] =
                        -alpha_eqn / (2.0 * (i - Nx));
                }
            }

            // ------------------------------------------------------------
            // v-pressure block
            // ------------------------------------------------------------
            if ((i > Nx + 2) && (i < 2 * Nx + 1) && (j > 2 * Nx))
            {
                if ((j < (i + Nx)) && (j == i + Nx - 1))
                {
                    A[i0][j0] = -(c_n_1 / gamma);
                }
                else if ((j > (i + Nx)) && (j == i + Nx + 1))
                {
                    A[i0][j0] = 1.0 / gamma;
                }
            }

            // ------------------------------------------------------------
            // Continuity equation block
            // ------------------------------------------------------------

            // u contribution to continuity
            // This block is removed in the special k = 0 matrix.
            if (!use_k0_matrix)
            {
                if ((i > 2 * Nx) && (j < Nx + 1))
                {
                    if (j == i - 2 * Nx)
                    {
                        A[i0][j0] = I * kx_current;
                    }
                    else if (j == i - 2 * Nx + 2)
                    {
                        A[i0][j0] = -I * kx_current;
                    }
                }
            }

            // v contribution to continuity
            if ((i > 2 * Nx)
                && (i < 3 * Nx - 1)
                && (j > Nx)
                && (j < 2 * Nx + 1))
            {
                if (j == i - Nx + 1)
                {
                    A[i0][j0] = 2.0 * (j - Nx - 1);
                }
            }
            else if ((i == 3 * Nx - 1)
                && (j > Nx)
                && (j < 2 * Nx + 1))
            {
                if (j == i - Nx + 1)
                {
                    A[i0][j0] = j - Nx - 1;
                }
            }
        }
    }

    // Extra conditions for k = 0 matrix
    if (use_k0_matrix)
    {
        // 1) Gauge condition: d0 = 0
        int gauge_row = Nx + 2;
        int col_d0 = 2 * Nx + 1;

        for (int j = 1; j <= Ntot; j++)
        {
            A[gauge_row - 1][j - 1] = Complex(0.0, 0.0);
        }

        A[gauge_row - 1][col_d0 - 1] = Complex(1.0, 0.0);
        b[gauge_row - 1] = Complex(0.0, 0.0);

        // 2) pressure BC: p_y(1) = 0
        int pbc_row = 3 * Nx;

        for (int j = 1; j <= Ntot; j++)
        {
            A[pbc_row - 1][j - 1] = Complex(0.0, 0.0);
        }

        b[pbc_row - 1] = Complex(0.0, 0.0);

        for (int n = 1; n <= Nx - 1; n++)
        {
            int col_dn = 2 * Nx + 1 + n;
            A[pbc_row - 1][col_dn - 1] = static_cast<double>(n * n);
        }
    }

    return solve_dense_linear_system(A, b);
}

UVPResult solve_uvp_chebyshev(const ChannelGrid& grid,
    const ComplexMatrix& q_n_u,
    const ComplexMatrix& q_n_v,
    double beta,
    double gamma)
{
    int Nx = grid.Nx;
    int Ny = grid.Nx;

    ComplexMatrix cheb_coeffs(3 * Nx, std::vector<Complex>(Nx, Complex(0.0, 0.0)));

    for (int ix = 0; ix < Nx; ix++)
    {
        bool use_k0_matrix = (ix == 0);

        ComplexVector m = solve_global_matrix_for_mode(
            grid,
            ix,
            q_n_u,
            q_n_v,
            beta,
            gamma,
            use_k0_matrix
        );

        for (int row = 0; row < 3 * Nx; row++)
        {
            cheb_coeffs[row][ix] = m[row];
        }
    }

    ComplexMatrix a_n = zero_complex_matrix_uvp(Ny, Nx);
    ComplexMatrix b_n = zero_complex_matrix_uvp(Ny, Nx);
    ComplexMatrix d_n = zero_complex_matrix_uvp(Ny, Nx);

    for (int row = 0; row < 3 * Nx; row++)
    {
        for (int ix = 0; ix < Nx; ix++)
        {
            if (row < Nx)
            {
                a_n[row][ix] = cheb_coeffs[row][ix];
            }
            else if ((row >= Nx) && (row < 2 * Nx))
            {
                b_n[row - Nx][ix] = cheb_coeffs[row][ix];
            }
            else
            {
                d_n[row - 2 * Nx][ix] = cheb_coeffs[row][ix];
            }
        }
    }

    ComplexMatrix u_hatk = zero_complex_matrix_uvp(Ny, Nx);
    ComplexMatrix v_hatk = zero_complex_matrix_uvp(Ny, Nx);
    ComplexMatrix P_hatk = zero_complex_matrix_uvp(Ny, Nx);

    for (int ix = 0; ix < Nx; ix++)
    {
        ComplexVector a_column(Ny);
        ComplexVector b_column(Ny);
        ComplexVector d_column(Ny);

        for (int n = 0; n < Ny; n++)
        {
            a_column[n] = a_n[n][ix];
            b_column[n] = b_n[n][ix];
            d_column[n] = d_n[n][ix];
        }

        ComplexVector u_column = cheb_eval_series(a_column, grid.y);
        ComplexVector v_column = cheb_eval_series(b_column, grid.y);
        ComplexVector P_column = cheb_eval_series(d_column, grid.y);

        for (int iy = 0; iy < Ny; iy++)
        {
            u_hatk[iy][ix] = u_column[iy];
            v_hatk[iy][ix] = v_column[iy];
            P_hatk[iy][ix] = P_column[iy];
        }
    }

    RealMatrix u = real_part_matrix_uvp(
        ifft_rows_complex_to_complex_fftw(u_hatk)
    );

    RealMatrix v = real_part_matrix_uvp(
        ifft_rows_complex_to_complex_fftw(v_hatk)
    );

    ComplexMatrix dPdx_hatk = zero_complex_matrix_uvp(Ny, Nx);
    Complex I(0.0, 1.0);

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int ix = 0; ix < Nx; ix++)
        {
            dPdx_hatk[iy][ix] = I * grid.k_x[ix] * P_hatk[iy][ix];
        }
    }

    RealMatrix dPdx = real_part_matrix_uvp(
        ifft_rows_complex_to_complex_fftw(dPdx_hatk)
    );

    ComplexMatrix dPdy_hatk = zero_complex_matrix_uvp(Ny, Nx);

    for (int ix = 0; ix < Nx; ix++)
    {
        ComplexVector P_column(Ny);

        for (int iy = 0; iy < Ny; iy++)
        {
            P_column[iy] = P_hatk[iy][ix];
        }

        ComplexVector a_P = a_cheb_coeff(P_column, grid.theta);
        ComplexVector b_P = b_cheb_coeff(a_P, Ny - 1);
        ComplexVector dPdy_column = cheb_eval_series(b_P, grid.y);

        for (int iy = 0; iy < Ny; iy++)
        {
            dPdy_hatk[iy][ix] = dPdy_column[iy];
        }
    }

    RealMatrix dPdy = real_part_matrix_uvp(
        ifft_rows_complex_to_complex_fftw(dPdy_hatk)
    );

    UVPResult result;

    result.u_hatk = u_hatk;
    result.v_hatk = v_hatk;
    result.P_hatk = P_hatk;

    result.u = u;
    result.v = v;

    result.dPdx = dPdx;
    result.dPdy = dPdy;

    return result;
}