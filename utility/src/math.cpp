#include "utility/math.h"
#include "utility/vectors.h"

namespace AIS4104::utility {
    // Help function for evaluating floating point numbers
    bool float_equals(const double a, const double b, const double epsilon = 1e-6)
    {
        return std::abs(a-b) < epsilon;
    }

    double cot(const double x)
    {
        return 1 / (std::sin(x) / std::cos(x)); // Use the known property that tan = sin/cos
    }

    // Equations from MR-preprint 2019 Appendix B chapter B.1.1
    Eigen::Vector3d euler_zyx_from_rotation_matrix(const Eigen::Matrix3d &r)
    {
        Eigen::Vector3d euler_zyx; // Init the vector for returning the rotations
        double alpha,beta,gamma;
        const double r31 = std::clamp(r(2,0), -1.0, 1.0); // Clamp to avoid domain issues from numerical drift

        if (!float_equals(r31, -1.0) && !float_equals(r31, 1.0)) {
            alpha = std::atan2(r(1, 0), r(0, 0));
            beta = std::atan2(-r31, std::hypot(r(0,0), r(1,0))); // Switched to using hypot instead of manually solving the root and square
            gamma = std::atan2(r(2,1), r(2,2));
        }
        else if (float_equals(r31, -1.0)) {
            alpha = 0.0;
            beta = std::numbers::pi / 2;
            gamma = std::atan2(r(0,1), r(1,1));
        }
        else if (float_equals(r31, 1.0)) {
            alpha = 0.0;
            beta = -std::numbers::pi / 2;
            gamma = -std::atan2(r(0,1), r(1,1));
        }

        euler_zyx << alpha, beta, gamma;

        return euler_zyx;
    }

    // Equation (3.30) page 74, MR pre-print 2019
    Eigen::Matrix3d skew_symmetric(const Eigen::Vector3d &v)
    {
        Eigen::Matrix3d skew;
        skew << 0.0, -v.z(), v.y(),
                v.z(), 0.0, -v.x(),
                -v.y(), v.x(), 0.0;
        return skew;
    }

    // Equation (3.30) page 74, MR pre-print 2019
    Eigen::Vector3d from_skew_symmetric(const Eigen::Matrix3d &m)
    {
            return { m(2,1), m(0,2), m(1,0) };
    }

    // Definition 3.20 on page 98, MR pre-print 2019
    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
    {
        Eigen::Matrix<double,6,6> Adj;
        Adj.setZero();
        Adj.topLeftCorner<3,3>() = r;
        Adj.bottomLeftCorner<3,3>() = skew_symmetric(p) * r;
        Adj.bottomRightCorner<3,3>() = r;
        return Adj;
    }

    // Definition 3.20 on page 98, MR pre-print 2019
    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf)
    {
        const Eigen::Matrix3d R = tf.topLeftCorner<3,3>();
        const Eigen::Vector3d p = tf.topRightCorner<3,1>();

        return adjoint_matrix(R, p);
    }

    // Definition 3.20 on page 98, MR pre-print 2019
    Eigen::VectorXd adjoint_map(const Eigen::VectorXd &twist, const Eigen::Matrix4d &tf)
    {
        return adjoint_matrix(tf) * twist;
    }

    // Equation 3.70 page 96, MR pre-print 2019
    Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
    {
        Eigen::VectorXd twist(6);
        twist << w, v;
        return twist;
    }

    // Unnamed equation page 101, MR pre-print 2019
    Eigen::VectorXd twist(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h, double angular_velocity)
    {
        // Angular part of the twist
        const Eigen::Vector3d w = angular_velocity * s;

        // Linear component
        const Eigen::Vector3d v = -angular_velocity * s.cross(q) + h * angular_velocity * s;
        return twist(w,v);
    }

    // Equation 3.71 page 96, MR pre-print 2019
    Eigen::Matrix4d twist_matrix(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
    {
        Eigen::Matrix4d twist_matrix = Eigen::Matrix4d::Zero();
        twist_matrix.topLeftCorner<3,3>() = skew_symmetric(w);
        twist_matrix.topRightCorner<3,1>() = v;
        return twist_matrix;
    }

    // Equation 3.71 page 96, MR pre-print 2019
    Eigen::Matrix4d twist_matrix(const Eigen::VectorXd &twist)
    {
        const Eigen::Vector3d w = twist.head<3>();
        const Eigen::Vector3d v = twist.tail<3>();
        return twist_matrix(w, v);
    }

    // Definition 3.24 page 102, MR pre-print 2019
    Eigen::VectorXd screw_axis(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
    {
        const Eigen::VectorXd Twist = twist(w, v);
        const double w_norm = w.norm();
        const double v_norm = v.norm();

        if (w_norm < 1e-6) {
            Eigen::Vector3d Screw = Twist / v_norm;
            return Screw;
        }

        Eigen::Vector3d Screw = Twist / w_norm;
        return Screw;
    }

    // Unnamed equation page 101, MR pre-print 2019
    Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h)
    {
        // Twist without angular velocity
        Eigen::Vector3d v = -s.cross(q) + h * s;
        Eigen::VectorXd axis(6);
        axis << s, v;
        return axis;
    }

    // Equation 3.51 page 82, MR pre-print 2019
    Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d &w, double theta)
    {
        const Eigen::Matrix3d w_hat = skew_symmetric(w);
        Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + std::sin(theta) * w_hat
                            + (1-std::cos(theta)) * w_hat * w_hat;
        return R;
    }

    // Proposition 3.25 page 103-104, MR pre-print 2019
    Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta)
    {
        // Test it we have a pure twist translation
        if (const double w_norm = w.norm(); w_norm < 1e-6) {
            Eigen::Vector3d p = v*theta;
            return transformation_matrix(p);
        }

        const auto R = matrix_exponential(w, theta);

        const Eigen::Matrix3d w_hat = skew_symmetric(w);
        const double c = std::cos(theta);
        const double s = std::sin(theta);

        const Eigen::Vector3d p = (Eigen::Matrix3d::Identity() * theta + (1-c)*w_hat
                                    + (theta-s)*w_hat*w_hat) * v;

        return transformation_matrix(R,p);
    }

    // Proposition 3.25 page 103-104, MR pre-print 2019
    Eigen::Matrix4d matrix_exponential(const Eigen::VectorXd &screw, double theta)
    {
        const Eigen::Vector3d w = screw.head<3>();
        const Eigen::Vector3d v = screw.tail<3>();
        return matrix_exponential(w,v,theta);
    }

    // Equations 3.58 - 3.61, pages 85 - 86, MR pre-print 2019
    std::pair<Eigen::Vector3d, double> matrix_logarithm(const Eigen::Matrix3d &r)
    {
        double theta;
        Eigen::Vector3d w;

        if (r.isApprox(Eigen::Matrix3d::Identity())) {
            theta = 0;
        }
        else {
            double trace_r = r.trace(); // Built-in eigen function for trace
            if (float_equals(trace_r, -1)) {
                theta = std::numbers::pi;
                w = (1 / sqrt(2 * (1 + r(2,2)))) * Eigen::Vector3d (r(0,2), r(1,2), 1 + r(2,2));
            }
            else {
                theta = std::acos(0.5 * (trace_r - 1));
                double w_n = 1.0 / 2.0 * std::sin(theta);

                double w_1 = w_n * (r(2, 1) - r(1, 2));
                double w_2 = w_n * (r(0, 2) - r(2, 0));
                double w_3 = w_n * (r(1, 0) - r(0, 1));

                w = Eigen::Vector3d(w_1, w_2, w_3);
            }
        }

        return {w, theta};
    }

    // Algorithm on page 104, MR pre-print 2019
    std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
    {
        // This is an idiotic solution, but a solution nonetheless
        Eigen::Matrix4d tf = transformation_matrix(r,p);
        return matrix_logarithm(tf);
    }

    // Algorithm on page 104, MR pre-print 2019
    std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix4d &tf)
    {
        const Eigen::Matrix3d R = tf.topLeftCorner<3,3>();
        const Eigen::Vector3d p = tf.topRightCorner<3,1>();

        Eigen::Vector3d w;
        Eigen::Vector3d v;
        double theta;

        if (R.isApprox(Eigen::Matrix3d::Identity())) {
            w = Eigen::Vector3d::Zero();
            v = p / p.norm();
            theta = p.norm();
        }
        else {
            std::pair<Eigen::Vector3d, double> m_log = matrix_logarithm(R);
            w = m_log.first;
            theta = m_log.second;
            const Eigen::Matrix3d skew_w = skew_symmetric(w);
            v = (((1/theta) * Eigen::Matrix3d::Identity()) - 0.5 * skew_w + ((1/theta) - 0.5* cot(theta/2)) * skew_w*skew_w) * p;
        }

        return {screw_axis(w,v), theta};
    }

    // Unnamed equation, page 72, MR pre-print 2019
    Eigen::Matrix3d rotate_x(double radians)
    {
        const double c = std::cos(radians);
        const double s = std::sin(radians);
        Eigen::Matrix3d rot_matrix;
        rot_matrix << 1.0, 0.0, 0.0,
                      0.0, c,-s,
                      0.0, s, c;
        return rot_matrix;
    }

    // Unnamed equation, page 72, MR pre-print 2019
    Eigen::Matrix3d rotate_y(double radians)
    {
        const double c = std::cos(radians);
        const double s = std::sin(radians);
        Eigen::Matrix3d rot_matrix;
        rot_matrix << c, 0.0, s,
                      0.0, 1.0, 0.0,
                     -s, 0.0, c;
        return rot_matrix;
    }

    // Unnamed equation, page 72, MR pre-print 2019
    Eigen::Matrix3d rotate_z(double radians)
    {
        const double c = std::cos(radians);
        const double s = std::sin(radians);
        Eigen::Matrix3d rot_matrix;
        rot_matrix << c,-s, 0.0,
                      s, c, 0.0,
                      0.0, 0.0, 1.0;
        return rot_matrix;
    }

    // Equations 3.17 page 67, MR pre-print 2019.
    Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x, const Eigen::Vector3d &y, const Eigen::Vector3d &z)
    {
        Eigen::Matrix3d rotation;
        rotation << x.normalized(), y.normalized(), z.normalized();
        return rotation;
    }

    // Euler rotation is just a way of creating a rotation matrix with a set order. Appendix B MR pre-print 2019
    Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e)
    {
        const double rot_x = e[2];
        const double rot_y = e[1];
        const double rot_z = e[0];

        Eigen::Matrix3d rot_matrix = rotate_z(rot_z) * rotate_y(rot_y) * rotate_x(rot_x);
        return rot_matrix;
    }

    // Equation 3.51 page 82, MR pre-print 2019
    // Rodrigues formula
    Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double radians)
    {
        const double c = std::cos(radians);
        const double s = std::sin(radians);

        const double n = axis.norm();

        const Eigen::Vector3d w = axis / n;
        Eigen::Matrix3d W = skew_symmetric(w);

        return Eigen::Matrix3d::Identity() + s * W + (1-c) * (W*W);
    }

    // Equation 3.62 page 87, MR pre-print 2019
    // Simple extraction of rotation from a homogenous transformation matrix
    Eigen::Matrix3d rotation_matrix(const Eigen::Matrix4d &tf)
    {
        return tf.topLeftCorner<3,3>();
    }

    // Equation 3.62 page 87, MR pre-print 2019
    // Creating a homogenous transformation matrix with only translation
    Eigen::Matrix4d transformation_matrix(const Eigen::Vector3d &p)
    {
        Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
        T.topRightCorner<3,1>() = p;
        return T;
    }

    // Equation 3.62 page 87, MR pre-print 2019
    // Creating a homogenous transformation matrix with only rotation
    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r)
    {
        Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
        T.topLeftCorner<3,3>() = r;
        return Eigen::Matrix4d::Identity();
    }

    // Equation 3.62 page 87, MR pre-print 2019
    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
    {
        Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
        T.topRightCorner<3,1>() = p;
        T.topLeftCorner<3,3>() = r;
        return T;
    }

}
