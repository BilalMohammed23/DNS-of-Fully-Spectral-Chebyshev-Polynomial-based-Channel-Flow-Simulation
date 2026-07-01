#include "chebyshev.h"

#include <cmath>

std::vector<double> a_cheb_coeff(const std::vector<double>& uhatk,
    const std::vector<double>& theta)
{
    int N = static_cast<int>(uhatk.size()) - 1;
    std::vector<double> a(N + 1, 0.0);

    for (int n = 0; n <= N; n++)
    {
        double sum = 0.0;

        for (int j = 0; j <= N; j++)
        {
            double c_j = 1.0;

            if (j == 0 || j == N)
            {
                c_j = 2.0;
            }

            sum += (1.0 / c_j) * uhatk[j] * std::cos(n * theta[j]);
        }

        double c_n = 1.0;

        if (n == 0 || n == N)
        {
            c_n = 2.0;
        }

        a[n] = (2.0 / (c_n * N)) * sum;
    }

    return a;
}

std::vector<double> b_cheb_coeff(const std::vector<double>& a,
    int N)
{
    std::vector<double> b(N + 1, 0.0);

    b[N] = 0.0;
    b[N - 1] = 2.0 * N * a[N];

    for (int k = N - 2; k >= 0; k--)
    {
        b[k] = b[k + 2] + 2.0 * (k + 1) * a[k + 1];
    }

    b[0] = b[0] / 2.0;

    return b;
}

std::vector<double> cheb_eval_series(const std::vector<double>& b,
    const std::vector<double>& x)
{
    int N = static_cast<int>(b.size()) - 1;
    std::vector<double> dudy(x.size(), 0.0);

    for (size_t i = 0; i < x.size(); i++)
    {
        double T0 = 1.0;

        if (N == 0)
        {
            dudy[i] = b[0] * T0;
            continue;
        }

        double T1 = x[i];
        dudy[i] = b[0] * T0 + b[1] * T1;

        for (int n = 1; n <= N - 1; n++)
        {
            double T2 = 2.0 * x[i] * T1 - T0;
            dudy[i] += b[n + 1] * T2;

            T0 = T1;
            T1 = T2;
        }
    }

    return dudy;
}

std::vector<Complex> a_cheb_coeff(const std::vector<Complex>& uhatk,
    const std::vector<double>& theta)
{
    int N = static_cast<int>(uhatk.size()) - 1;
    std::vector<Complex> a(N + 1, Complex(0.0, 0.0));

    for (int n = 0; n <= N; n++)
    {
        Complex sum(0.0, 0.0);

        for (int j = 0; j <= N; j++)
        {
            double c_j = 1.0;

            if (j == 0 || j == N)
            {
                c_j = 2.0;
            }

            sum += (1.0 / c_j) * uhatk[j] * std::cos(n * theta[j]);
        }

        double c_n = 1.0;

        if (n == 0 || n == N)
        {
            c_n = 2.0;
        }

        a[n] = (2.0 / (c_n * N)) * sum;
    }

    return a;
}

std::vector<Complex> b_cheb_coeff(const std::vector<Complex>& a,
    int N)
{
    std::vector<Complex> b(N + 1, Complex(0.0, 0.0));

    b[N] = Complex(0.0, 0.0);
    b[N - 1] = 2.0 * static_cast<double>(N) * a[N];

    for (int k = N - 2; k >= 0; k--)
    {
        b[k] = b[k + 2] + 2.0 * static_cast<double>(k + 1) * a[k + 1];
    }

    b[0] = b[0] / 2.0;

    return b;
}

std::vector<Complex> cheb_eval_series(const std::vector<Complex>& b,
    const std::vector<double>& x)
{
    int N = static_cast<int>(b.size()) - 1;
    std::vector<Complex> dudy(x.size(), Complex(0.0, 0.0));

    for (size_t i = 0; i < x.size(); i++)
    {
        double T0 = 1.0;

        if (N == 0)
        {
            dudy[i] = b[0] * T0;
            continue;
        }

        double T1 = x[i];
        dudy[i] = b[0] * T0 + b[1] * T1;

        for (int n = 1; n <= N - 1; n++)
        {
            double T2 = 2.0 * x[i] * T1 - T0;
            dudy[i] += b[n + 1] * T2;

            T0 = T1;
            T1 = T2;
        }
    }

    return dudy;
}