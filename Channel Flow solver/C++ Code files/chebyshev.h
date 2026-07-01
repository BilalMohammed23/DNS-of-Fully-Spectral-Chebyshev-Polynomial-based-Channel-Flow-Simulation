#pragma once

#include <vector>
#include <complex>

using Complex = std::complex<double>;

std::vector<double> a_cheb_coeff(const std::vector<double>& uhatk,
    const std::vector<double>& theta);

std::vector<double> b_cheb_coeff(const std::vector<double>& a,
    int N);

std::vector<double> cheb_eval_series(const std::vector<double>& b,
    const std::vector<double>& x);


// Complex-valued versions for Fourier-Chebyshev solver
std::vector<Complex> a_cheb_coeff(const std::vector<Complex>& uhatk,
    const std::vector<double>& theta);

std::vector<Complex> b_cheb_coeff(const std::vector<Complex>& a,
    int N);

std::vector<Complex> cheb_eval_series(const std::vector<Complex>& b,
    const std::vector<double>& x);