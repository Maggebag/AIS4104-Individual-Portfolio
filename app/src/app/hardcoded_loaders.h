#ifndef AIS4104_HARDCODED_LOADERS_H
#define AIS4104_HARDCODED_LOADERS_H

#include "app/examples/tracikkinematicssolver.h"

#include "app/implementations/screwskinematicssolver.h"

#include <utility/math.h>
#include <utility/vectors.h>

#include <memory>

namespace AIS4104 {

inline std::shared_ptr<ScrewsKinematicsSolver> hardcoded_ur3e_screw_solver()
{
    double w1 = 0.13105;
    double w2 = 0.0921;
    double l1 = 0.24355;
    double l2 = 0.2132;
    double h1 = 0.15185;
    double h2 = 0.08535;

    Eigen::Matrix4d m = utility::transformation_matrix(utility::rotate_y(-90.0 * utility::deg_to_rad) * utility::rotate_x(-90.0 * utility::deg_to_rad) * utility::rotate_z(-90.0 * utility::deg_to_rad),
        Eigen::Vector3d{l1 + l2, w1 + w2, h1 - h2});

    Simulation::JointLimits limits
    {
        utility::to_eigen_vectord(std::vector<double>{180.0, 180.0, 180.0, 360.0, 360.0, 360.0}) * utility::deg_to_rad,
        utility::to_eigen_vectord(std::vector<double>{90.0, 90.0, 90.0, 180.0, 180.0, 180.0}) * utility::deg_to_rad,
        utility::to_eigen_vectord(std::vector<double>{-360.0, -360.0, -360.0, -360.0, -360.0, -360.0}) * utility::deg_to_rad,
        utility::to_eigen_vectord(std::vector<double>{360.0, 360.0, 360.0, 360.0, 360.0, 360.0}) * utility::deg_to_rad
    };

    return std::make_shared<ScrewsKinematicsSolver>(
        m,
        std::vector<Eigen::VectorXd>{
            utility::screw_axis({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 0.0),
            utility::screw_axis({0.0, 0.0, h1}, {0.0, 1.0, 0.0}, 0.0),
            utility::screw_axis({l1, 0.0, h1}, {0.0, 1.0, 0.0}, 0.0),
            utility::screw_axis({l1 + l2, 0.0, h1}, {0.0, 1.0, 0.0}, 0.0),
            utility::screw_axis({l1 + l2, w1, 0.0}, {0.0, 0.0, -1.0}, 0.0),
            utility::screw_axis({l1 + l2, 0, h1 - h2}, {0.0, 1.0, 0.0}, 0.0)
        }, limits
    );
}

inline std::shared_ptr<TracIkKinematicsSolver> hardcoded_ur3e_tracik_solver()
{
    double h1 = 0.15185;
    double h2 = 0.08535;
    double w1 = 0.13105;
    double w2 = 0.0921;
    double l1 = 0.24355;
    double l2 = 0.2132;

    auto c = std::make_shared<KDL::Chain>();
    c->addSegment(KDL::Segment(KDL::Joint(KDL::Joint::JointType::Fixed), KDL::Frame(KDL::Vector(0.0, 0.0, h1))));
    c->addSegment(KDL::Segment(KDL::Joint(KDL::Joint::JointType::RotZ), KDL::Frame(KDL::Vector(0.0, w1, 0.0))));
    c->addSegment(KDL::Segment(KDL::Joint(KDL::Joint::JointType::RotY), KDL::Frame(KDL::Vector(l1, 0.0, 0.0))));
    c->addSegment(KDL::Segment(KDL::Joint(KDL::Joint::JointType::RotY), KDL::Frame(KDL::Vector(l2, 0.0, 0.0))));
    c->addSegment(KDL::Segment(KDL::Joint(KDL::Joint::JointType::RotY), KDL::Frame(KDL::Vector(0.0, 0.0, 0.0))));
    c->addSegment(KDL::Segment(KDL::Joint(KDL::Joint::JointType::RotY), KDL::Frame(KDL::Rotation::RotY(EIGEN_PI), KDL::Vector(0.0, 0.0, -h2))));
    c->addSegment(KDL::Segment(KDL::Joint(KDL::Joint::JointType::RotZ), KDL::Frame(KDL::Rotation::RotX(-0.5 * EIGEN_PI), KDL::Vector(0.0, w2, 0.0))));

    Simulation::JointLimits limits
    {
        utility::to_eigen_vectord(std::vector<double>{180.0, 180.0, 180.0, 360.0, 360.0, 360.0}) * utility::deg_to_rad,
        utility::to_eigen_vectord(std::vector<double>{90.0, 90.0, 90.0, 180.0, 180.0, 180.0}) * utility::deg_to_rad,
        utility::to_eigen_vectord(std::vector<double>{-360.0, -360.0, -360.0, -360.0, -360.0, -360.0}) * utility::deg_to_rad,
        utility::to_eigen_vectord(std::vector<double>{360.0, 360.0, 360.0, 360.0, 360.0, 360.0}) * utility::deg_to_rad
    };

    return std::make_shared<TracIkKinematicsSolver>(c, limits);
}

//TASK: Kinematic modeling of the KUKA KR 6 r900 sixx using screws.
inline std::shared_ptr<ScrewsKinematicsSolver> hardcoded_kr6r_screw_solver()
{
    return hardcoded_ur3e_screw_solver();
}

//TASK: Kinematic modeling of the KUKA KR 6 r900 sixx using Trac IK and KDL.
inline std::shared_ptr<TracIkKinematicsSolver> hardcoded_kr6r_tracik_solver()
{
    return hardcoded_ur3e_tracik_solver();
}

}

#endif
