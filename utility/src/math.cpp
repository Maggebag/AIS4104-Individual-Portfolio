#include "utility/math.h"
#include "utility/vectors.h"

namespace AIS4104::utility {

Eigen::Vector3d euler_zyx_from_rotation_matrix(const Eigen::Matrix3d &r)
{
    return Eigen::Vector3d::Zero();
}

Eigen::Matrix3d skew_symmetric(const Eigen::Vector3d &v)
{
    return Eigen::Matrix3d::Identity();
}

Eigen::Vector3d from_skew_symmetric(const Eigen::Matrix3d &m)
{
    return Eigen::Vector3d::Zero();
}

Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
{
    return Eigen::MatrixXd::Identity(6, 6);
}

Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf)
{
    return Eigen::MatrixXd::Identity(6, 6);
}

Eigen::VectorXd adjoint_map(const Eigen::VectorXd &twist, const Eigen::Matrix4d &tf)
{
    return Eigen::VectorXd::Zero(6);
}

Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
{
    return Eigen::VectorXd::Zero(6);
}

Eigen::VectorXd twist(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h, double angular_velocity)
{
    return Eigen::VectorXd::Zero(6);
}

Eigen::Matrix4d twist_matrix(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
{
    return Eigen::Matrix4d::Identity();
}

Eigen::Matrix4d twist_matrix(const Eigen::VectorXd &twist)
{
    return Eigen::Matrix4d::Identity();
}

Eigen::VectorXd screw_axis(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
{
    return Eigen::VectorXd::Zero(6);
}

Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h)
{
    return Eigen::VectorXd::Zero(6);
}

Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d &w, double theta)
{
    return Eigen::Matrix3d::Identity();
}

Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta)
{
    return Eigen::Matrix4d::Identity();
}

Eigen::Matrix4d matrix_exponential(const Eigen::VectorXd &screw, double theta)
{
    return Eigen::Matrix4d::Identity();
}

std::pair<Eigen::Vector3d, double> matrix_logarithm(const Eigen::Matrix3d &r)
{
    return std::make_pair(Eigen::Vector3d::Zero(), 0);
}

std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
{
    return std::make_pair(Eigen::VectorXd::Zero(6), 0);
}

std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix4d &tf)
{
    return std::make_pair(Eigen::VectorXd::Zero(6), 0);
}

Eigen::Matrix3d rotate_x(double radians)
{
    return Eigen::Matrix3d::Identity();
}

Eigen::Matrix3d rotate_y(double radians)
{
    return Eigen::Matrix3d::Identity();
}

Eigen::Matrix3d rotate_z(double radians)
{
    return Eigen::Matrix3d::Identity();
}

Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x, const Eigen::Vector3d &y, const Eigen::Vector3d &z)
{
    return Eigen::Matrix3d::Identity();
}

Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e)
{
    return Eigen::Matrix3d::Identity();
}

Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double radians)
{
    return Eigen::Matrix3d::Identity();
}

Eigen::Matrix3d rotation_matrix(const Eigen::Matrix4d &tf)
{
    return Eigen::Matrix3d::Identity();
}

Eigen::Matrix4d transformation_matrix(const Eigen::Vector3d &p)
{
    return Eigen::Matrix4d::Identity();
}

Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r)
{
    return Eigen::Matrix4d::Identity();
}

Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
{
    return Eigen::Matrix4d::Identity();
}

}
