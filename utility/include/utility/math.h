#ifndef AIS4104_UTILITY_MATH_H
#define AIS4104_UTILITY_MATH_H

#include <Eigen/Dense>
#include <numbers>

namespace AIS4104::utility {

constexpr double rad_to_deg = 57.2957795;
constexpr double deg_to_rad = 0.0174532925;

//TASK: Implement the corresponding function definitions in utility/src/math.cpp

Eigen::Vector3d euler_zyx_from_rotation_matrix(const Eigen::Matrix3d &r);

Eigen::Matrix3d skew_symmetric(const Eigen::Vector3d &v);
Eigen::Vector3d from_skew_symmetric(const Eigen::Matrix3d &m);

Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p);
Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf);
Eigen::VectorXd adjoint_map(const Eigen::VectorXd &twist, const Eigen::Matrix4d &tf);

Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v);
Eigen::VectorXd twist(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h, double angular_velocity);

Eigen::Matrix4d twist_matrix(const Eigen::Vector3d &w, const Eigen::Vector3d &v);
Eigen::Matrix4d twist_matrix(const Eigen::VectorXd &twist);

Eigen::VectorXd screw_axis(const Eigen::Vector3d &w, const Eigen::Vector3d &v);
Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h);

Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d &w, double theta);
Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta);
Eigen::Matrix4d matrix_exponential(const Eigen::VectorXd &screw, double theta);

std::pair<Eigen::Vector3d, double> matrix_logarithm(const Eigen::Matrix3d &r);
std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix3d &r, const Eigen::Vector3d &p);
std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix4d &tf);

Eigen::Matrix3d rotate_x(double radians);
Eigen::Matrix3d rotate_y(double radians);
Eigen::Matrix3d rotate_z(double radians);
Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x, const Eigen::Vector3d &y, const Eigen::Vector3d &z);
Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e);
Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double radians);

Eigen::Matrix3d rotation_matrix(const Eigen::Matrix4d &tf);
Eigen::Matrix4d transformation_matrix(const Eigen::Vector3d &p);
Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r);
Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p);

}

#endif
