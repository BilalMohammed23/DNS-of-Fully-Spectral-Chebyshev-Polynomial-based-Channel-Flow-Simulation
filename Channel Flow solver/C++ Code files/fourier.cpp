#include "fourier.h"

#include <fftw3.h>
#include <stdexcept>

ComplexMatrix fft_rows_real_to_complex_fftw(const RealMatrix& input)
{
    if (input.empty() || input[0].empty())
    {
        throw std::runtime_error("Input matrix is empty.");
    }

    int Ny = static_cast<int>(input.size());
    int Nx = static_cast<int>(input[0].size());

    ComplexMatrix output(Ny, std::vector<Complex>(Nx, Complex(0.0, 0.0)));

    fftw_complex* in = static_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * Nx)
        );

    fftw_complex* out = static_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * Nx)
        );

    if (!in || !out)
    {
        fftw_free(in);
        fftw_free(out);
        throw std::runtime_error("FFTW memory allocation failed.");
    }

    fftw_plan plan = fftw_plan_dft_1d(
        Nx,
        in,
        out,
        FFTW_FORWARD,
        FFTW_ESTIMATE
    );

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int ix = 0; ix < Nx; ix++)
        {
            in[ix][0] = input[iy][ix];
            in[ix][1] = 0.0;
        }

        fftw_execute(plan);

        for (int k = 0; k < Nx; k++)
        {
            output[iy][k] = Complex(out[k][0] / static_cast<double>(Nx),
                out[k][1] / static_cast<double>(Nx));
        }
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return output;
}

RealMatrix ifft_rows_complex_to_real_fftw(const ComplexMatrix& input)
{
    if (input.empty() || input[0].empty())
    {
        throw std::runtime_error("Input matrix is empty.");
    }

    int Ny = static_cast<int>(input.size());
    int Nx = static_cast<int>(input[0].size());

    RealMatrix output(Ny, std::vector<double>(Nx, 0.0));

    fftw_complex* in = static_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * Nx)
        );

    fftw_complex* out = static_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * Nx)
        );

    if (!in || !out)
    {
        fftw_free(in);
        fftw_free(out);
        throw std::runtime_error("FFTW memory allocation failed.");
    }

    fftw_plan plan = fftw_plan_dft_1d(
        Nx,
        in,
        out,
        FFTW_BACKWARD,
        FFTW_ESTIMATE
    );

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int k = 0; k < Nx; k++)
        {
            in[k][0] = input[iy][k].real();
            in[k][1] = input[iy][k].imag();
        }

        fftw_execute(plan);

        for (int ix = 0; ix < Nx; ix++)
        {
            output[iy][ix] = out[ix][0];
        }
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return output;
}

ComplexMatrix fft_rows_complex_to_complex_fftw(const ComplexMatrix& input)
{
    if (input.empty() || input[0].empty())
    {
        throw std::runtime_error("Input matrix is empty.");
    }

    int Ny = static_cast<int>(input.size());
    int Nx = static_cast<int>(input[0].size());

    ComplexMatrix output(Ny, std::vector<Complex>(Nx, Complex(0.0, 0.0)));

    fftw_complex* in = static_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * Nx)
        );

    fftw_complex* out = static_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * Nx)
        );

    if (!in || !out)
    {
        fftw_free(in);
        fftw_free(out);
        throw std::runtime_error("FFTW memory allocation failed.");
    }

    fftw_plan plan = fftw_plan_dft_1d(
        Nx,
        in,
        out,
        FFTW_FORWARD,
        FFTW_ESTIMATE
    );

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int ix = 0; ix < Nx; ix++)
        {
            in[ix][0] = input[iy][ix].real();
            in[ix][1] = input[iy][ix].imag();
        }

        fftw_execute(plan);

        for (int k = 0; k < Nx; k++)
        {
            output[iy][k] = Complex(out[k][0] / static_cast<double>(Nx),
                out[k][1] / static_cast<double>(Nx));
        }
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return output;
}

ComplexMatrix ifft_rows_complex_to_complex_fftw(const ComplexMatrix& input)
{
    if (input.empty() || input[0].empty())
    {
        throw std::runtime_error("Input matrix is empty.");
    }

    int Ny = static_cast<int>(input.size());
    int Nx = static_cast<int>(input[0].size());

    ComplexMatrix output(Ny, std::vector<Complex>(Nx, Complex(0.0, 0.0)));

    fftw_complex* in = static_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * Nx)
        );

    fftw_complex* out = static_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * Nx)
        );

    if (!in || !out)
    {
        fftw_free(in);
        fftw_free(out);
        throw std::runtime_error("FFTW memory allocation failed.");
    }

    fftw_plan plan = fftw_plan_dft_1d(
        Nx,
        in,
        out,
        FFTW_BACKWARD,
        FFTW_ESTIMATE
    );

    for (int iy = 0; iy < Ny; iy++)
    {
        for (int k = 0; k < Nx; k++)
        {
            in[k][0] = input[iy][k].real();
            in[k][1] = input[iy][k].imag();
        }

        fftw_execute(plan);

        for (int ix = 0; ix < Nx; ix++)
        {
            output[iy][ix] = Complex(out[ix][0], out[ix][1]);
        }
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return output;
}