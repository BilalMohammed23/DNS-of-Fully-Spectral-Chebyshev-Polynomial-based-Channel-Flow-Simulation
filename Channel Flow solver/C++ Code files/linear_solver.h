#pragma once

#include <vector>
#include <complex>

using Complex = std::complex<double>;
using ComplexVector = std::vector<Complex>;
using ComplexDenseMatrix = std::vector<std::vector<Complex>>;

ComplexVector solve_dense_linear_system(ComplexDenseMatrix A,
    ComplexVector b);