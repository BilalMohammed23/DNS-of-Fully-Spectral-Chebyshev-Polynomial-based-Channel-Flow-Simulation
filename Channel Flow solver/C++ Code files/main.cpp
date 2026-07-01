#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

#include "channel_solver.h"
#include "output.h"

static double max_abs_matrix_main(const RealMatrix& A)
{
    double max_value = 0.0;

    for (size_t iy = 0; iy < A.size(); iy++)
    {
        for (size_t ix = 0; ix < A[iy].size(); ix++)
        {
            max_value = std::max(max_value, std::abs(A[iy][ix]));
        }
    }

    return max_value;
}

int main()
{
    double T = 10;

    ChannelSolverResult result = run_channel_solver(T);

    std::cout << std::endl;
    std::cout << std::setprecision(12);
    std::cout << "Full solver completed." << std::endl;
    std::cout << "Final time = " << result.final_time << std::endl;
    std::cout << "Step count = " << result.step_count << std::endl;
    std::cout << "Final max divergence = " << result.max_divergence_final << std::endl;
    std::cout << "Final max |u| = " << max_abs_matrix_main(result.u_final) << std::endl;
    std::cout << "Final max |v| = " << max_abs_matrix_main(result.v_final) << std::endl;

    write_vector_csv("x_grid.csv", result.grid.x);
    write_vector_csv("y_grid.csv", result.grid.y);

    write_matrix_csv("u_final.csv", result.u_final);
    write_matrix_csv("v_final.csv", result.v_final);

    write_centerline_profiles("centerline_profiles.csv",
        result.grid,
        result.u_final,
        result.v_final);

    std::cout << std::endl;
    std::cout << "CSV files written successfully." << std::endl;

    return 0;
}