#include "linear_solver.h"

#include <stdexcept>
#include <cmath>
#include <algorithm>

ComplexVector solve_dense_linear_system(ComplexDenseMatrix A,
    ComplexVector b)
{
    int n = static_cast<int>(A.size());

    if (n == 0 || static_cast<int>(b.size()) != n)
    {
        throw std::runtime_error("Invalid linear system size.");
    }

    for (int i = 0; i < n; i++)
    {
        if (static_cast<int>(A[i].size()) != n)
        {
            throw std::runtime_error("Matrix must be square.");
        }
    }

    // Gaussian elimination with partial pivoting
    for (int k = 0; k < n; k++)
    {
        int pivot_row = k;
        double pivot_abs = std::abs(A[k][k]);

        for (int i = k + 1; i < n; i++)
        {
            double candidate_abs = std::abs(A[i][k]);

            if (candidate_abs > pivot_abs)
            {
                pivot_abs = candidate_abs;
                pivot_row = i;
            }
        }

        if (pivot_abs < 1.0e-14)
        {
            throw std::runtime_error("Singular or nearly singular matrix in solve_dense_linear_system.");
        }

        if (pivot_row != k)
        {
            std::swap(A[k], A[pivot_row]);
            std::swap(b[k], b[pivot_row]);
        }

        for (int i = k + 1; i < n; i++)
        {
            Complex factor = A[i][k] / A[k][k];

            for (int j = k; j < n; j++)
            {
                A[i][j] -= factor * A[k][j];
            }

            b[i] -= factor * b[k];
        }
    }

    ComplexVector x(n, Complex(0.0, 0.0));

    // Back substitution
    for (int i = n - 1; i >= 0; i--)
    {
        Complex sum = b[i];

        for (int j = i + 1; j < n; j++)
        {
            sum -= A[i][j] * x[j];
        }

        x[i] = sum / A[i][i];
    }

    return x;
}