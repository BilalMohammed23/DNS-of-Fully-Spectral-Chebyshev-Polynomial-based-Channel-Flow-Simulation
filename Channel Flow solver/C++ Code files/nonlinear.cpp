#include "nonlinear.h"
#include "chebyshev.h"

#include <complex>
#include <stdexcept>

static ComplexMatrix zero_complex_matrix(int Ny, int Nx)
{
    return ComplexMatrix(Ny, std::vector<Complex>(Nx, Complex(0.0, 0.0)));
}

static RealMatrix real_part_matrix(const ComplexMatrix& input)
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

static ComplexMatrix derivative_x_hat(const ComplexMatrix& fhat,
    const std::vector<double>& k_x)
{
    int Ny = static_cast<int>(fhat.size());
    int Nx = static_cast<int>(fhat[0].size());

    ComplexMatrix dfdx_hat = zero_complex_matrix(Ny, Nx);

    Complex I(0.0, 1.0);

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int k = 0; k < Nx; k++)
        {
            dfdx_hat[iy][k] = I * k_x[k] * fhat[iy][k];
        }
    }

    return dfdx_hat;
}

static ComplexMatrix derivative_y_hat(const ComplexMatrix& fhat,
    const std::vector<double>& theta,
    const std::vector<double>& y)
{
    int Ny = static_cast<int>(fhat.size());
    int Nx = static_cast<int>(fhat[0].size());

    ComplexMatrix dfdy_hat = zero_complex_matrix(Ny, Nx);

    for (int k = 0; k < Nx; k++)
    {
        std::vector<Complex> column(Ny);

        for (int iy = 0; iy < Ny; iy++)
        {
            column[iy] = fhat[iy][k];
        }

        std::vector<Complex> a = a_cheb_coeff(column, theta);
        std::vector<Complex> b = b_cheb_coeff(a, Ny - 1);
        std::vector<Complex> dfdy_col = cheb_eval_series(b, y);

        for (int iy = 0; iy < Ny; iy++)
        {
            dfdy_hat[iy][k] = dfdy_col[iy];
        }
    }

    return dfdy_hat;
}

static ComplexMatrix pad_to_three_halves(const ComplexMatrix& input)
{
    int Ny = static_cast<int>(input.size());
    int Nx = static_cast<int>(input[0].size());
    int Nx_pssp = 3 * Nx / 2;

    ComplexMatrix output = zero_complex_matrix(Ny, Nx_pssp);

    // positive modes: MATLAB columns 1:Nx/2
    for (int iy = 0; iy < Ny; iy++)
    {
        for (int k = 0; k < Nx / 2; k++)
        {
            output[iy][k] = input[iy][k];
        }
    }

    // negative modes:
    // MATLAB: output(:, Nx+1:end) = input(:, Nx/2+1:end)
    // C++ zero-based: output[:, Nx : Nx_pssp-1] = input[:, Nx/2 : Nx-1]
    for (int iy = 0; iy < Ny; iy++)
    {
        for (int k = Nx / 2; k < Nx; k++)
        {
            int kp = k + Nx / 2;
            output[iy][kp] = input[iy][k];
        }
    }

    return output;
}

static ComplexMatrix truncate_from_three_halves(const ComplexMatrix& input_padded,
    int Nx_original)
{
    int Ny = static_cast<int>(input_padded.size());
    int Nx_pssp = static_cast<int>(input_padded[0].size());

    ComplexMatrix output = zero_complex_matrix(Ny, Nx_original);

    // positive modes
    for (int iy = 0; iy < Ny; iy++)
    {
        for (int k = 0; k < Nx_original / 2; k++)
        {
            output[iy][k] = input_padded[iy][k];
        }
    }

    // negative modes:
    // MATLAB: output(:, Nx/2+1:end) = input_padded(:, Nx+1:end)
    // C++ zero-based: output[:, Nx/2 : Nx-1] = input_padded[:, Nx_original : Nx_pssp-1]
    for (int iy = 0; iy < Ny; iy++)
    {
        for (int k = Nx_original / 2; k < Nx_original; k++)
        {
            int kp = k + Nx_original / 2;
            if (kp >= Nx_pssp)
            {
                throw std::runtime_error("Invalid truncation index.");
            }

            output[iy][k] = input_padded[iy][kp];
        }
    }

    return output;
}

static ComplexMatrix multiply_fields(const ComplexMatrix& a,
    const ComplexMatrix& b)
{
    int Ny = static_cast<int>(a.size());
    int Nx = static_cast<int>(a[0].size());

    ComplexMatrix output = zero_complex_matrix(Ny, Nx);

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int ix = 0; ix < Nx; ix++)
        {
            output[iy][ix] = a[iy][ix] * b[iy][ix];
        }
    }

    return output;
}

static RealMatrix add_real_matrices(const RealMatrix& a,
    const RealMatrix& b)
{
    int Ny = static_cast<int>(a.size());
    int Nx = static_cast<int>(a[0].size());

    RealMatrix output(Ny, std::vector<double>(Nx, 0.0));

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int ix = 0; ix < Nx; ix++)
        {
            output[iy][ix] = a[iy][ix] + b[iy][ix];
        }
    }

    return output;
}

NonlinearResult compute_nonlinear_terms(const ChannelGrid& grid,
    const ComplexMatrix& uhatk,
    const ComplexMatrix& vhatk)
{
    int Ny = static_cast<int>(uhatk.size());
    int Nx = static_cast<int>(uhatk[0].size());

    if (Nx != grid.Nx || Ny != grid.Nx)
    {
        throw std::runtime_error("Size mismatch in compute_nonlinear_terms.");
    }

    // Fourier and Chebyshev derivatives in spectral space
    ComplexMatrix dudx_hatk = derivative_x_hat(uhatk, grid.k_x);
    ComplexMatrix dudy_hatk = derivative_y_hat(uhatk, grid.theta, grid.y);

    ComplexMatrix dvdx_hatk = derivative_x_hat(vhatk, grid.k_x);
    ComplexMatrix dvdy_hatk = derivative_y_hat(vhatk, grid.theta, grid.y);

    // 3/2 padding
    ComplexMatrix uhatk_pssp = pad_to_three_halves(uhatk);
    ComplexMatrix vhatk_pssp = pad_to_three_halves(vhatk);

    ComplexMatrix dudx_hatk_pssp = pad_to_three_halves(dudx_hatk);
    ComplexMatrix dudy_hatk_pssp = pad_to_three_halves(dudy_hatk);

    ComplexMatrix dvdx_hatk_pssp = pad_to_three_halves(dvdx_hatk);
    ComplexMatrix dvdy_hatk_pssp = pad_to_three_halves(dvdy_hatk);

    // Back to physical pseudo-spectral grid
    ComplexMatrix u_pssp = ifft_rows_complex_to_complex_fftw(uhatk_pssp);
    ComplexMatrix v_pssp = ifft_rows_complex_to_complex_fftw(vhatk_pssp);

    ComplexMatrix dudx_pssp = ifft_rows_complex_to_complex_fftw(dudx_hatk_pssp);
    ComplexMatrix dudy_pssp = ifft_rows_complex_to_complex_fftw(dudy_hatk_pssp);

    ComplexMatrix dvdx_pssp = ifft_rows_complex_to_complex_fftw(dvdx_hatk_pssp);
    ComplexMatrix dvdy_pssp = ifft_rows_complex_to_complex_fftw(dvdy_hatk_pssp);

    // Nonlinear products on padded physical grid
    ComplexMatrix H1_1_pssp = multiply_fields(u_pssp, dudx_pssp);
    ComplexMatrix H1_2_pssp = multiply_fields(v_pssp, dudy_pssp);

    ComplexMatrix H2_1_pssp = multiply_fields(u_pssp, dvdx_pssp);
    ComplexMatrix H2_2_pssp = multiply_fields(v_pssp, dvdy_pssp);

    // FFT of nonlinear products
    ComplexMatrix H1_1_pssp_f = fft_rows_complex_to_complex_fftw(H1_1_pssp);
    ComplexMatrix H1_2_pssp_f = fft_rows_complex_to_complex_fftw(H1_2_pssp);

    ComplexMatrix H2_1_pssp_f = fft_rows_complex_to_complex_fftw(H2_1_pssp);
    ComplexMatrix H2_2_pssp_f = fft_rows_complex_to_complex_fftw(H2_2_pssp);

    // Truncate back to Nx spectral grid
    ComplexMatrix postprod_H1_1 = truncate_from_three_halves(H1_1_pssp_f, Nx);
    ComplexMatrix postprod_H1_2 = truncate_from_three_halves(H1_2_pssp_f, Nx);

    ComplexMatrix postprod_H2_1 = truncate_from_three_halves(H2_1_pssp_f, Nx);
    ComplexMatrix postprod_H2_2 = truncate_from_three_halves(H2_2_pssp_f, Nx);

    // Back to original physical grid
    RealMatrix H1_1 = real_part_matrix(ifft_rows_complex_to_complex_fftw(postprod_H1_1));
    RealMatrix H1_2 = real_part_matrix(ifft_rows_complex_to_complex_fftw(postprod_H1_2));

    RealMatrix H2_1 = real_part_matrix(ifft_rows_complex_to_complex_fftw(postprod_H2_1));
    RealMatrix H2_2 = real_part_matrix(ifft_rows_complex_to_complex_fftw(postprod_H2_2));

    NonlinearResult result;

    result.H1 = add_real_matrices(H1_1, H1_2);
    result.H2 = add_real_matrices(H2_1, H2_2);

    result.dudx_hatk = dudx_hatk;
    result.dvdy_hatk = dvdy_hatk;

    return result;
}