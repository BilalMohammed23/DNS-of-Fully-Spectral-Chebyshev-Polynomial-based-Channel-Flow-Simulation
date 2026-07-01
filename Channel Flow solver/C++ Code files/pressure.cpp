#include "pressure.h"
#include "chebyshev.h"
#include "linear_solver.h"

#include <stdexcept>
#include <cmath>

static ComplexMatrix zero_complex_matrix_pressure(int Ny, int Nx)
{
    return ComplexMatrix(Ny, std::vector<Complex>(Nx, Complex(0.0, 0.0)));
}

static RealMatrix zero_real_matrix_pressure(int Ny, int Nx)
{
    return RealMatrix(Ny, std::vector<double>(Nx, 0.0));
}

static RealMatrix real_part_matrix_pressure(const ComplexMatrix& input)
{
    int Ny = static_cast<int>(input.size());
    int Nx = static_cast<int>(input[0].size());

    RealMatrix output = zero_real_matrix_pressure(Ny, Nx);

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int ix = 0; ix < Nx; ix++)
        {
            output[iy][ix] = input[iy][ix].real();
        }
    }

    return output;
}

static ComplexVector solve_tridiag_dn_pressure(int Nx,
    const ComplexVector& q_n,
    double alpha,
    Complex ghat,
    Complex fhat)
{
    int half = Nx / 2;

    ComplexDenseMatrix A_even(half, ComplexVector(half, Complex(0.0, 0.0)));
    ComplexDenseMatrix A_odd(half, ComplexVector(half, Complex(0.0, 0.0)));

    ComplexVector b_even(half, Complex(0.0, 0.0));
    ComplexVector b_odd(half, Complex(0.0, 0.0));

    double c_n = 1.0;
    double c_n_1 = 1.0;

    // even modes: n = 0, 2, 4, ..., Nx-2
    for (int n = 0; n <= Nx - 2; n += 2)
    {
        int row = n / 2;

        double c_n_2 = 1.0;
        if (n - 2 == 0)
        {
            c_n_2 = 2.0;
        }

        if (alpha == 0.0)
        {
            if (n == 0)
            {
                b_even[row] = ghat - fhat;
            }
            else if (n == 2)
            {
                b_even[row] = Complex(0.0, 0.0);
            }
            else if (n == Nx - 2)
            {
                b_even[row] =
                    ((c_n_1 * q_n[n] - c_n_1 * c_n_2 * q_n[n - 2]) / (2.0 * (n - 1)))
                    + ((c_n * q_n[n]) / (2.0 * (n + 1)));
            }
            else
            {
                b_even[row] =
                    ((c_n_1 * q_n[n] - c_n_1 * c_n_2 * q_n[n - 2]) / (2.0 * (n - 1)))
                    + ((c_n * q_n[n] - q_n[n + 2]) / (2.0 * (n + 1)));
            }
        }
        else
        {
            if (n == 0)
            {
                b_even[row] = ghat - fhat;
            }
            else if (n == Nx - 2)
            {
                b_even[row] =
                    ((c_n_1 * q_n[n] - c_n_1 * c_n_2 * q_n[n - 2]) / (2.0 * (n - 1)))
                    + ((c_n * q_n[n]) / (2.0 * (n + 1)));
            }
            else
            {
                b_even[row] =
                    ((c_n_1 * q_n[n] - c_n_1 * c_n_2 * q_n[n - 2]) / (2.0 * (n - 1)))
                    + ((c_n * q_n[n] - q_n[n + 2]) / (2.0 * (n + 1)));
            }
        }

        for (int col = 0; col < half; col++)
        {
            if (alpha == 0.0)
            {
                if (row == 0)
                {
                    double mode = static_cast<double>(2 * col);
                    A_even[row][col] = 2.0 * mode * mode;
                }
                else if (row == 1 && col == 0)
                {
                    A_even[row][col] = 1.0;
                }

                if ((row == col) && (row > 1))
                {
                    A_even[row][col] =
                        -(2.0 * n
                            + ((alpha * c_n_1) / (2.0 * (n - 1)))
                            + ((alpha * c_n) / (2.0 * (n + 1))));
                }
                else if ((row > 1) && (col < row) && (col == row - 1))
                {
                    A_even[row][col] = (alpha * c_n_2 * c_n_1) / (2.0 * (n - 1));
                }
                else if ((row > 1) && (col > row) && (col == row + 1))
                {
                    A_even[row][col] = alpha / (2.0 * (n + 1));
                }
            }
            else
            {
                if (row == 0)
                {
                    double mode = static_cast<double>(2 * col);
                    A_even[row][col] = 2.0 * mode * mode;
                }

                if ((row == col) && (row > 0))
                {
                    A_even[row][col] =
                        -(2.0 * n
                            + ((alpha * c_n_1) / (2.0 * (n - 1)))
                            + ((alpha * c_n) / (2.0 * (n + 1))));
                }
                else if ((row > 0) && (col < row) && (col == row - 1))
                {
                    A_even[row][col] = (alpha * c_n_2 * c_n_1) / (2.0 * (n - 1));
                }
                else if ((row > 0) && (col > row) && (col == row + 1))
                {
                    A_even[row][col] = alpha / (2.0 * (n + 1));
                }
            }
        }
    }

    // odd modes: n = 1, 3, 5, ..., Nx-1
    for (int n = 1; n <= Nx - 1; n += 2)
    {
        int row = (n - 1) / 2;

        double c_n_2 = 1.0;

        if (n == 1)
        {
            b_odd[row] = ghat + fhat;
        }
        else if (n == Nx - 1)
        {
            b_odd[row] =
                ((c_n_1 * q_n[n] - c_n_1 * c_n_2 * q_n[n - 2]) / (2.0 * (n - 1)))
                + ((c_n * q_n[n]) / (2.0 * (n + 1)));
        }
        else
        {
            b_odd[row] =
                ((c_n_1 * q_n[n] - c_n_1 * c_n_2 * q_n[n - 2]) / (2.0 * (n - 1)))
                + ((c_n * q_n[n] - q_n[n + 2]) / (2.0 * (n + 1)));
        }

        for (int col = 0; col < half; col++)
        {
            if (row == 0)
            {
                double mode = static_cast<double>(2 * col + 1);
                A_odd[row][col] = 2.0 * mode * mode;
            }

            if ((row == col) && (row > 0))
            {
                A_odd[row][col] =
                    -(2.0 * n
                        + ((alpha * c_n_1) / (2.0 * (n - 1)))
                        + ((alpha * c_n) / (2.0 * (n + 1))));
            }
            else if ((row > 0) && (col < row) && (col == row - 1))
            {
                A_odd[row][col] = (alpha * c_n_2 * c_n_1) / (2.0 * (n - 1));
            }
            else if ((row > 0) && (col > row) && (col == row + 1))
            {
                A_odd[row][col] = alpha / (2.0 * (n + 1));
            }
        }
    }

    ComplexVector a_even = solve_dense_linear_system(A_even, b_even);
    ComplexVector a_odd = solve_dense_linear_system(A_odd, b_odd);

    ComplexVector a_n(Nx, Complex(0.0, 0.0));

    for (int n = 0; n <= Nx - 2; n += 2)
    {
        a_n[n] = a_even[n / 2];
    }

    for (int n = 1; n <= Nx - 1; n += 2)
    {
        a_n[n] = a_odd[(n - 1) / 2];
    }

    return a_n;
}

PressureResult compute_pressure_0th_timestep(const ChannelGrid& grid,
    const RealMatrix& H1,
    const RealMatrix& H2)
{
    int Nx = grid.Nx;
    int Ny = grid.Nx;

    Complex ghat(0.0, 0.0);
    Complex fhat(0.0, 0.0);
    Complex I(0.0, 1.0);

    // dH1/dx
    ComplexMatrix H1_hatk = fft_rows_real_to_complex_fftw(H1);
    ComplexMatrix dH1dx_hatk = zero_complex_matrix_pressure(Ny, Nx);

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int k = 0; k < Nx; k++)
        {
            dH1dx_hatk[iy][k] = I * grid.k_x[k] * H1_hatk[iy][k];
        }
    }

    RealMatrix dH1dx = real_part_matrix_pressure(
        ifft_rows_complex_to_complex_fftw(dH1dx_hatk)
    );

    // dH2/dy in physical x-columns
    RealMatrix dH2dy = zero_real_matrix_pressure(Ny, Nx);

    for (int ix = 0; ix < Nx; ix++)
    {
        std::vector<double> column(Ny);

        for (int iy = 0; iy < Ny; iy++)
        {
            column[iy] = H2[iy][ix];
        }

        std::vector<double> a = a_cheb_coeff(column, grid.theta);
        std::vector<double> b = b_cheb_coeff(a, Ny - 1);
        std::vector<double> derivative_column = cheb_eval_series(b, grid.y);

        for (int iy = 0; iy < Ny; iy++)
        {
            dH2dy[iy][ix] = derivative_column[iy];
        }
    }

    // q = -(dH1dx + dH2dy)
    RealMatrix q = zero_real_matrix_pressure(Ny, Nx);

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int ix = 0; ix < Nx; ix++)
        {
            q[iy][ix] = -(dH1dx[iy][ix] + dH2dy[iy][ix]);
        }
    }

    ComplexMatrix q_hatk = fft_rows_real_to_complex_fftw(q);

    // Solve pressure Poisson equation mode-by-mode
    ComplexMatrix P_hatk = zero_complex_matrix_pressure(Ny, Nx);

    for (int k = 0; k < Nx; k++)
    {
        ComplexVector q_column(Ny);

        for (int iy = 0; iy < Ny; iy++)
        {
            q_column[iy] = q_hatk[iy][k];
        }

        ComplexVector q_n = a_cheb_coeff(q_column, grid.theta);

        double alpha = grid.k_x[k] * grid.k_x[k];

        ComplexVector d_n = solve_tridiag_dn_pressure(Nx, q_n, alpha, ghat, fhat);

        ComplexVector P_column = cheb_eval_series(d_n, grid.y);

        for (int iy = 0; iy < Ny; iy++)
        {
            P_hatk[iy][k] = P_column[iy];
        }
    }

    ComplexMatrix P_complex = ifft_rows_complex_to_complex_fftw(P_hatk);
    RealMatrix P = real_part_matrix_pressure(P_complex);

    // dP/dx
    ComplexMatrix dPdx_hatk = zero_complex_matrix_pressure(Ny, Nx);

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int k = 0; k < Nx; k++)
        {
            dPdx_hatk[iy][k] = I * grid.k_x[k] * P_hatk[iy][k];
        }
    }

    RealMatrix dPdx = real_part_matrix_pressure(
        ifft_rows_complex_to_complex_fftw(dPdx_hatk)
    );

    // dP/dy
    ComplexMatrix dPdy_hatk = zero_complex_matrix_pressure(Ny, Nx);

    for (int k = 0; k < Nx; k++)
    {
        ComplexVector P_column(Ny);

        for (int iy = 0; iy < Ny; iy++)
        {
            P_column[iy] = P_hatk[iy][k];
        }

        ComplexVector a = a_cheb_coeff(P_column, grid.theta);
        ComplexVector b = b_cheb_coeff(a, Ny - 1);
        ComplexVector derivative_column = cheb_eval_series(b, grid.y);

        for (int iy = 0; iy < Ny; iy++)
        {
            dPdy_hatk[iy][k] = derivative_column[iy];
        }
    }

    RealMatrix dPdy = real_part_matrix_pressure(
        ifft_rows_complex_to_complex_fftw(dPdy_hatk)
    );

    PressureResult result;
    result.P = P;
    result.dPdx = dPdx;
    result.dPdy = dPdy;

    return result;
}