#pragma once

#include <vector>
#include <complex>

using Complex = std::complex<double>;
using RealMatrix = std::vector<std::vector<double>>;
using ComplexMatrix = std::vector<std::vector<Complex>>;

ComplexMatrix fft_rows_real_to_complex_fftw(const RealMatrix& input);

RealMatrix ifft_rows_complex_to_real_fftw(const ComplexMatrix& input);

ComplexMatrix fft_rows_complex_to_complex_fftw(const ComplexMatrix& input);

ComplexMatrix ifft_rows_complex_to_complex_fftw(const ComplexMatrix& input);