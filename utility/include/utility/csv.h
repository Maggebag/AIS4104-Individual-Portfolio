#ifndef AIS4104_UTILITY_CSV_H
#define AIS4104_UTILITY_CSV_H

#include "utility/vectors.h"

#include <time.h>
#include <fstream>
#include <filesystem>

#include <Eigen/src/Core/IO.h>

namespace AIS4104::utility {

inline std::filesystem::path date_time_stamped_filename(const std::string &name)
{
    auto t = std::time(nullptr);
    std::ostringstream oss;
    oss << std::put_time(std::localtime(&t), "%Y-%m-%d_%H-%M-%S") << "_" << name;
    return oss.str();
}

template<typename F> requires utility::is_floating_point<F>
void eigen_vector_to_csv_file(const std::vector<Eigen::Matrix<F, Eigen::Dynamic, 1>> &vectors, const std::filesystem::path &path)
{
    std::ofstream csv_file(path);
    Eigen::IOFormat csv_format(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ");
    for(auto &vector : vectors)
        csv_file << vector.transpose().format(csv_format) << std::endl;
}

template<typename T> requires std::is_fundamental<T>::value
void std_vector_to_csv_file(const std::vector<T> &vector, const std::filesystem::path &path)
{
    std::ofstream csv_file(path);
    Eigen::IOFormat csv_format(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ");
    for(auto i = 0u; i < vector.size(); i++)
        csv_file << vector[i] << std::endl;
}

}

#endif
