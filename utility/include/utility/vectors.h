#ifndef AIS4104_UTILITY_VECTORS_H
#define AIS4104_UTILITY_VECTORS_H

#include <vector>

#include <kdl/jntarray.hpp>

namespace AIS4104::utility {

std::vector<float> to_std_vectorf(const Eigen::VectorXf &vector);
std::vector<float> to_std_vectorf(const Eigen::VectorXd &vector);

std::vector<double> to_std_vectord(const Eigen::VectorXd &vector);
std::vector<double> to_std_vectord(const Eigen::VectorXf &vector);

Eigen::VectorXd to_eigen_vectord(const std::vector<float> &vector);
Eigen::VectorXd to_eigen_vectord(const std::vector<double> &vector);

Eigen::VectorXf to_eigen_vectorf(const std::vector<float> &vector);
Eigen::VectorXf to_eigen_vectorf(const std::vector<double> &vector);

Eigen::VectorXd to_eigen_vectord(const Eigen::VectorXf &vector);
Eigen::VectorXf to_eigen_vectorf(const Eigen::VectorXd &vector);

KDL::JntArray to_kdl_array(const Eigen::VectorXd &vector);

Eigen::Matrix4d to_eigen_matrix4d(const KDL::Frame &f);
Eigen::Vector3d to_eigen_vector3d(const KDL::Vector &v);
Eigen::Matrix3d to_eigen_matrix3d(const KDL::Rotation &r);

KDL::Frame to_kdl_frame(const Eigen::Matrix4d &tf);
KDL::Vector to_kdl_vector(const Eigen::Vector3d &vector);
KDL::Rotation to_kdl_rotation(const Eigen::Matrix3d &rotation);

template<typename F>
concept is_floating_point = std::floating_point<std::remove_reference_t<F>>;

template<typename F1, typename F2, typename F = F1> requires is_floating_point<F1> && is_floating_point<F2> && is_floating_point<F>
bool is_approx_equal(F1 a, F2 b, F eps = std::numeric_limits<F>::epsilon())
{
    return std::fabs(a - b) < eps;
}

template<typename F1, typename F2, typename F = F1> requires is_floating_point<F1> && is_floating_point<F2> && is_floating_point<F>
bool is_approx_equal(const std::vector<F1> &first, const std::vector<F2> &second, F eps = std::numeric_limits<F>::epsilon())
{
    if(first.size() != second.size())
        return false;
    for(auto i = 0u; i < first.size(); i++)
        if(!is_approx_equal(first[i], second[i], eps))
            return false;
    return true;
}

template<typename F1, typename F2, typename F = F1> requires is_floating_point<F1> && is_floating_point<F2> && is_floating_point<F>
bool is_approx_equal(const std::vector<F1> &first, const Eigen::Matrix<F2, Eigen::Dynamic, 1> &second, F eps = std::numeric_limits<F>::epsilon())
{
    if(first.size() != second.size())
        return false;
    for(auto i = 0u; i < first.size(); i++)
        if(!is_approx_equal(first[i], second[i], eps))
            return false;
    return true;
}

template<typename F1, typename F2, typename F = F1> requires is_floating_point<F1> && is_floating_point<F2> && is_floating_point<F>
bool is_approx_equal(const Eigen::Matrix<F1, Eigen::Dynamic, 1> &first, const std::vector<F2> &second, F eps = std::numeric_limits<F>::epsilon())
{
    return is_approx_equal<F2, F1, F>(second, first, eps);
}

template<typename F1, typename F2, typename F = F1> requires is_floating_point<F1> && is_floating_point<F2> && is_floating_point<F>
bool is_approx_equal(const Eigen::Matrix<F1, Eigen::Dynamic, 1> &first, const Eigen::Matrix<F2, Eigen::Dynamic, 1> &second, F eps = std::numeric_limits<F>::epsilon())
{
    if(first.size() != second.size())
        return false;
    for(auto i = 0u; i < first.size(); i++)
        if(!is_approx_equal(first[i], second[i], eps))
            return false;
    return true;
}

template<typename F> requires is_floating_point<F>
bool is_approx_equal(const Eigen::Matrix<F, Eigen::Dynamic, 1> &first, const Eigen::Matrix<F, Eigen::Dynamic, 1> &second, F eps = std::numeric_limits<F>::epsilon())
{
    if(first.size() != second.size())
        return false;
    for(auto i = 0u; i < first.size(); i++)
        if(!is_approx_equal(first[i], second[i], eps))
            return false;
    return true;
}

template<typename F, int I> requires is_floating_point<F>
bool is_approx_equal(const Eigen::Matrix<F, I, 1> &first, const Eigen::Matrix<F, I, 1> &second, F eps = std::numeric_limits<F>::epsilon())
{
    for(auto i = 0u; i < I; i++)
        if(!is_approx_equal(first[i], second[i], eps))
            return false;
    return true;
}

template<typename F, int I, int J> requires is_floating_point<F>
bool is_approx_equal(const Eigen::Matrix<F, I, J> &first, const Eigen::Matrix<F, I, J> &second, F eps = std::numeric_limits<F>::epsilon())
{
    for(auto i = 0u; i < I; i++)
        for(auto j = 0u; j < J; j++)
            if(!is_approx_equal(first(i, j), second(i, j), eps))
                return false;
    return true;
}

inline bool is_approx_equal(const KDL::Vector &first, const KDL::Vector &second, double eps = std::numeric_limits<double>::epsilon())
{
    return is_approx_equal(first.data[0], second.data[0], eps)
        && is_approx_equal(first.data[1], second.data[1], eps)
        && is_approx_equal(first.data[2], second.data[2], eps);
}

template<typename F> requires is_floating_point<F>
bool is_approx_equal(const KDL::Vector &first, const Eigen::Matrix<F, 3, 1> &second, F eps = std::numeric_limits<F>::epsilon())
{
    return is_approx_equal(first.data[0], second[0], eps)
        && is_approx_equal(first.data[1], second[1], eps)
        && is_approx_equal(first.data[2], second[2], eps);
}

template<typename F> requires is_floating_point<F>
bool is_approx_equal(const Eigen::Matrix<F, 3, 1> &first, const KDL::Vector &second, F eps = std::numeric_limits<F>::epsilon())
{
    return is_approx_equal(second, first);
}

inline bool is_approx_equal(const KDL::Rotation &first, const KDL::Rotation &second, double eps = std::numeric_limits<double>::epsilon())
{
    for(auto i = 0u; i < 9; i++)
        if(!is_approx_equal(first.data[i], second.data[i], eps))
            return false;
    return true;
}

template<typename F> requires is_floating_point<F>
bool is_approx_equal(const KDL::Rotation &first, const Eigen::Matrix<F, 3, 3> &second, F eps = std::numeric_limits<F>::epsilon())
{
    for(auto i = 0u; i < 3; i++)
        for(auto j = 0u; j < 3; j++)
            if(!is_approx_equal(first(i, j), second(i, j), eps))
                return false;
    return true;
}

template<typename F> requires is_floating_point<F>
bool is_approx_equal(const Eigen::Matrix<F, 3, 3> &first, const KDL::Rotation &second, F eps = std::numeric_limits<F>::epsilon())
{
    return is_approx_equal(second, first, eps);
}

template<typename F> requires is_floating_point<F>
bool is_approx_equal(const Eigen::Matrix<F, 4, 4> &first, const KDL::Frame &second, F eps = std::numeric_limits<F>::epsilon())
{
    Eigen::Matrix3d r = first.block(0, 0, 3, 3);
    Eigen::Vector3d p = first.block(0, 3, 3, 1);
    return is_approx_equal(r, second.M, eps)
        && is_approx_equal(p, second.p, eps);
}

template<typename F> requires is_floating_point<F>
bool is_approx_equal(const KDL::Frame &first, const Eigen::Matrix<F, 4, 4> &second, F eps = std::numeric_limits<F>::epsilon())
{
    return is_approx_equal(second, first, eps);
}

inline bool is_approx_equal(const KDL::Frame &first, const KDL::Frame &second, double eps = std::numeric_limits<double>::epsilon())
{
    return is_approx_equal(first.p, second.p, eps)
        && is_approx_equal(first.M, second.M, eps);
}

}

#endif
