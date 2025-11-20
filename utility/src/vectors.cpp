#include "utility/vectors.h"

#include "utility/math.h"

namespace AIS4104::utility {

std::vector<float> to_std_vectorf(const Eigen::VectorXf &vector)
{
    std::vector<float> ret(vector.size());
    std::memcpy(ret.data(), vector.data(), vector.size() * sizeof(float));
    return ret;
}

std::vector<float> to_std_vectorf(const Eigen::VectorXd &vector)
{
    std::vector<float> ret(vector.size());
    for(auto i = 0u; i < vector.size(); i++)
        ret[i] = static_cast<float>(vector[i]);
    return ret;
}

std::vector<double> to_std_vectord(const Eigen::VectorXd &vector)
{
    std::vector<double> ret(vector.size());
    std::memcpy(ret.data(), vector.data(), vector.size() * sizeof(double));
    return ret;
}

std::vector<double> to_std_vectord(const Eigen::VectorXf &vector)
{
    std::vector<double> ret(vector.size());
    for(auto i = 0u; i < vector.size(); i++)
        ret[i] = vector[i];
    return ret;
}

Eigen::VectorXd to_eigen_vectord(const std::vector<float> &vector)
{
    Eigen::VectorXd ret(vector.size());
    for(auto i = 0u; i < vector.size(); i++)
        ret[i] = vector[i];
    return ret;
}

Eigen::VectorXd to_eigen_vectord(const std::vector<double> &vector)
{
    Eigen::VectorXd result(vector.size());
    std::memcpy(result.data(), vector.data(), vector.size() * sizeof(double));
    return result;
}

Eigen::VectorXf to_eigen_vectorf(const std::vector<float> &vector)
{
    Eigen::VectorXf ret(vector.size());
    std::memcpy(ret.data(), vector.data(), vector.size() * sizeof(float));
    return ret;
}

Eigen::VectorXf to_eigen_vectorf(const std::vector<double> &vector)
{
    Eigen::VectorXf ret(vector.size());
    for(auto i = 0u; i < vector.size(); i++)
        ret[i] = static_cast<float>(vector[i]);
    return ret;
}

Eigen::VectorXd to_eigen_vectord(const Eigen::VectorXf &vector)
{
    Eigen::VectorXd ret(vector.size());
    for(auto i = 0u; i < vector.size(); i++)
        ret[i] = static_cast<double>(vector[i]);
    return ret;
}

Eigen::VectorXf to_eigen_vectorf(const Eigen::VectorXd &vector)
{
    Eigen::VectorXf ret(vector.size());
    for(auto i = 0u; i < vector.size(); i++)
        ret[i] = static_cast<float>(vector[i]);
    return ret;
}

KDL::JntArray to_kdl_array(const Eigen::VectorXd &vector)
{
    KDL::JntArray ret(static_cast<uint32_t>(vector.size()));
    ret.data = vector;
    return ret;
}

KDL::Frame to_kdl_frame(const Eigen::Matrix4d &tf)
{
    return KDL::Frame(to_kdl_rotation(tf.block<3, 3>(0, 0)), to_kdl_vector(tf.block<3, 1>(0, 3)));
}

KDL::Vector to_kdl_vector(const Eigen::Vector3d &vector)
{
    KDL::Vector ret;
    std::memcpy(ret.data, vector.data(), vector.size() * sizeof(double));
    return ret;
}

KDL::Rotation to_kdl_rotation(const Eigen::Matrix3d &rotation)
{
    return KDL::Rotation(
        rotation(0, 0), rotation(0, 1), rotation(0, 2)
      , rotation(1, 0), rotation(1, 1), rotation(1, 2)
      , rotation(2, 0), rotation(2, 1), rotation(2, 2)
    );
}

Eigen::Matrix4d to_eigen_matrix4d(const KDL::Frame &f)
{
    return transformation_matrix(to_eigen_matrix3d(f.M), to_eigen_vector3d(f.p));
}

Eigen::Vector3d to_eigen_vector3d(const KDL::Vector &v)
{
    return {v.x(), v.y(), v.z()};
}

Eigen::Matrix3d to_eigen_matrix3d(const KDL::Rotation &r)
{
    Eigen::Matrix3d ret;
    ret <<
        r.data[0], r.data[1], r.data[2],
        r.data[3], r.data[4], r.data[5],
        r.data[6], r.data[7], r.data[8];
    return ret;
}

}
