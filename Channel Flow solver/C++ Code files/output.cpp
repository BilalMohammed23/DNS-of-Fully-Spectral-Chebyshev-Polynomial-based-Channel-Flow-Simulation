#include "output.h"

#include <fstream>
#include <stdexcept>

void write_vector_csv(const std::string& filename,
    const std::vector<double>& data)
{
    std::ofstream file(filename);

    if (!file)
    {
        throw std::runtime_error("Could not open file: " + filename);
    }

    for (size_t i = 0; i < data.size(); i++)
    {
        file << data[i];

        if (i + 1 < data.size())
        {
            file << ",";
        }
    }

    file << "\n";
}

void write_matrix_csv(const std::string& filename,
    const RealMatrix& data)
{
    std::ofstream file(filename);

    if (!file)
    {
        throw std::runtime_error("Could not open file: " + filename);
    }

    for (size_t i = 0; i < data.size(); i++)
    {
        for (size_t j = 0; j < data[i].size(); j++)
        {
            file << data[i][j];

            if (j + 1 < data[i].size())
            {
                file << ",";
            }
        }

        file << "\n";
    }
}

void write_centerline_profiles(const std::string& filename,
    const ChannelGrid& grid,
    const RealMatrix& u,
    const RealMatrix& v)
{
    std::ofstream file(filename);

    if (!file)
    {
        throw std::runtime_error("Could not open file: " + filename);
    }

    int Nx = grid.Nx;
    int ix_center = Nx / 2;

    file << "y,u_center,v_center\n";

    for (int iy = 0; iy < Nx; iy++)
    {
        file << grid.y[iy] << ","
            << u[iy][ix_center] << ","
            << v[iy][ix_center] << "\n";
    }
}