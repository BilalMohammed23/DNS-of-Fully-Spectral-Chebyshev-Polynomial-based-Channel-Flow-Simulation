#include "viscous.h"
#include "chebyshev.h"

#include <stdexcept>

static ComplexMatrix zero_complex_matrix_viscous(int Ny, int Nx)
{
    return ComplexMatrix(Ny, std::vector<Complex>(Nx, Complex(0.0, 0.0)));
}

static RealMatrix real_part_matrix_viscous(const ComplexMatrix& input)
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

static ComplexMatrix second_derivative_x_hat(const ComplexMatrix& fhat,
    const std::vector<double>& k_x)
{
    int Ny = static_cast<int>(fhat.size());
    int Nx = static_cast<int>(fhat[0].size());

    ComplexMatrix d2fdx2_hat = zero_complex_matrix_viscous(Ny, Nx);

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int k = 0; k < Nx; k++)
        {
            d2fdx2_hat[iy][k] = -(k_x[k] * k_x[k]) * fhat[iy][k];
        }
    }

    return d2fdx2_hat;
}

static ComplexMatrix second_derivative_y_hat(const ComplexMatrix& fhat,
    const std::vector<double>& theta,
    const std::vector<double>& y)
{
    int Ny = static_cast<int>(fhat.size());
    int Nx = static_cast<int>(fhat[0].size());

    ComplexMatrix d2fdy2_hat = zero_complex_matrix_viscous(Ny, Nx);

    for (int k = 0; k < Nx; k++)
    {
        std::vector<Complex> column(Ny);

        for (int iy = 0; iy < Ny; iy++)
        {
            column[iy] = fhat[iy][k];
        }

        std::vector<Complex> a = a_cheb_coeff(column, theta);
        std::vector<Complex> b = b_cheb_coeff(a, Ny - 1);
        std::vector<Complex> c = b_cheb_coeff(b, Ny - 1);
        std::vector<Complex> d2fdy2_col = cheb_eval_series(c, y);

        for (int iy = 0; iy < Ny; iy++)
        {
            d2fdy2_hat[iy][k] = d2fdy2_col[iy];
        }
    }

    return d2fdy2_hat;
}

ViscousResult compute_viscous_terms(const ChannelGrid& grid,
    const ComplexMatrix& uhatk,
    const ComplexMatrix& vhatk)
{
    int Ny = static_cast<int>(uhatk.size());
    int Nx = static_cast<int>(uhatk[0].size());

    if (Nx != grid.Nx || Ny != grid.Nx)
    {
        throw std::runtime_error("Size mismatch in compute_viscous_terms.");
    }

    ComplexMatrix d2udx2_hatk = second_derivative_x_hat(uhatk, grid.k_x);
    ComplexMatrix d2vdx2_hatk = second_derivative_x_hat(vhatk, grid.k_x);

    ComplexMatrix d2udy2_hatk = second_derivative_y_hat(uhatk, grid.theta, grid.y);
    ComplexMatrix d2vdy2_hatk = second_derivative_y_hat(vhatk, grid.theta, grid.y);

    ViscousResult result;

    result.d2udx2 = real_part_matrix_viscous(
        ifft_rows_complex_to_complex_fftw(d2udx2_hatk)
    );

    result.d2vdx2 = real_part_matrix_viscous(
        ifft_rows_complex_to_complex_fftw(d2vdx2_hatk)
    );

    result.d2udy2 = real_part_matrix_viscous(
        ifft_rows_complex_to_complex_fftw(d2udy2_hatk)
    );

    result.d2vdy2 = real_part_matrix_viscous(
        ifft_rows_complex_to_complex_fftw(d2vdy2_hatk)
    );

    return result;
}