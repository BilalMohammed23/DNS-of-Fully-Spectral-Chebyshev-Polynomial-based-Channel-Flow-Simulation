#pragma once

#include <string>
#include "grid.h"
#include "fourier.h"

void write_vector_csv(const std::string& filename,
    const std::vector<double>& data);

void write_matrix_csv(const std::string& filename,
    const RealMatrix& data);

void write_centerline_profiles(const std::string& filename,
    const ChannelGrid& grid,
    const RealMatrix& u,
    const RealMatrix& v);