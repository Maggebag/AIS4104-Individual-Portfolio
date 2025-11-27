#include "app/robotwrapper.h"

#include <utility/math.h>
#include <utility/vectors.h>
#include <iostream>

using namespace AIS4104;

RobotWrapper::RobotWrapper(std::shared_ptr<threepp::Robot> robot, std::shared_ptr<Simulation::KinematicsSolver> solver)
    : m_tool_transform(Eigen::Matrix4d::Identity())
    , m_robot(std::move(robot))
    , m_solver(std::move(solver))
{
}

threepp::Robot& RobotWrapper::threepp_robot()
{
    return *m_robot;
}

std::shared_ptr<threepp::Robot> RobotWrapper::threepp_robot_ptr()
{
    return m_robot;
}

const threepp::Robot& RobotWrapper::threepp_robot() const
{
    return *m_robot;
}

uint8_t RobotWrapper::joint_count() const
{
    return m_solver->joint_count();
}

//Completed-TASK: Implement the function to calculate the joint positions for the desired tool pose
// a) Use m_tool_transform to calculate the flange pose required by m_solver.ik_solve()
// b) Use the m_solver.ik_solve() overload with the solution selector lambda to choose the most desirable IK solution.
Eigen::VectorXd RobotWrapper::ik_solve_pose(const Eigen::Matrix4d &eef_pose, const Eigen::VectorXd &j0) const
{
    Eigen::Matrix4d T_tool_flange = m_tool_transform.inverse();
    Eigen::Matrix4d desired_flange_pose = eef_pose * T_tool_flange;
    return ik_solve_flange_pose(desired_flange_pose, j0);
}

//Completed-TASK: Implement the function to calculate the joint positions for the desired flange pose
// a) Use m_tool_transform to calculate the flange pose required by m_solver.ik_solve()
// b) Use the m_solver.ik_solve() overload with the solution selector lambda to choose the most desirable IK solution.
Eigen::VectorXd RobotWrapper::ik_solve_flange_pose(const Eigen::Matrix4d &flange_pose, const Eigen::VectorXd &j0) const
{
    // Lambda solution selector, pick solution with the smallest angle change in joints
    auto selector = [&](const std::vector<Eigen::VectorXd> &candidates)->uint32_t {
        uint32_t best_idx = 0;
        double best_cost = 0.0;

        auto lower_limits = m_solver->joint_limits().lower_position;
        auto upper_limits = m_solver->joint_limits().upper_position;

        double limit_margin = 0.1;
        double weight_distance = 1.0;
        double weight_limit = 10.0;

        for (uint32_t i = 0; i < candidates.size(); i++) {
            const auto &q = candidates[i];
            double cost = 0.0;

            for (int j = 0; j < q.size(); j++) {
                cost += weight_distance * std::abs(q(j) - j0(j));

                double dist_to_min = std::abs(q(j) - lower_limits(j));
                double dist_to_max = std::abs(upper_limits(j) - q(j));

                double d = std::min(dist_to_min, dist_to_max);

                if (d < limit_margin) { cost += weight_limit * (limit_margin - d); }
            }
            if (cost < best_cost) {
                best_cost = cost;
                best_idx = i;
            }
        }
        return best_idx;
    };

    return m_solver->ik_solve(flange_pose, j0, selector);
}

Eigen::Matrix4d RobotWrapper::tool_transform() const
{
    return m_tool_transform;
}

void RobotWrapper::set_tool_transform(Eigen::Matrix4d transform)
{
    m_tool_transform = std::move(transform);
}

//TASK: Calculate the pose of the end effector using forward kinematics;
// Relevant variables are m_solver and m_tool_transform.
Eigen::Matrix4d RobotWrapper::current_pose() const
{
    Eigen::Matrix4d T = current_flange_pose();
    return T * m_tool_transform;
}

//TASK: Calculate the position of the end effector using forward kinematics.
// Relevant variables are m_solver and m_tool_transform (or possibly another function of RobotWrapper?).
Eigen::Vector3d RobotWrapper::current_position() const
{
    Eigen::Matrix4d T = current_pose();
    return T.topRightCorner<3,1>();
}

//TASK: Calculate the orientation of the end effector using forward kinematics and m_solver (or rely on another function of RobotWrapper?).
Eigen::Vector3d RobotWrapper::current_orientation_zyx() const
{
    Eigen::Matrix4d T = current_pose();
    Eigen::Matrix3d R = utility::rotation_matrix(T);
    return utility::euler_zyx_from_rotation_matrix(R);
}

//TASK: Calculate the pose of the end effector using forward kinematics and m_solver.
Eigen::Matrix4d RobotWrapper::current_flange_pose() const
{
    Eigen::VectorXd q = joint_positions();
    return m_solver->fk_solve(q);
}

//TASK: Based on the flange pose, return its linear position.
Eigen::Vector3d RobotWrapper::current_flange_position() const
{
    Eigen::Matrix4d T = current_pose();
    return T.topRightCorner<3, 1>();
}

//TASK: Based on the flange pose, return its orientation in the Euler ZYX representation.
Eigen::Vector3d RobotWrapper::current_flange_orientation_zyx() const
{
    Eigen::Matrix4d T = current_pose();
    Eigen::Matrix3d R = utility::rotation_matrix(T);
    return utility::euler_zyx_from_rotation_matrix(R);
}

const Simulation::JointLimits& RobotWrapper::joint_limits() const
{
    return m_solver->joint_limits();
}

Eigen::VectorXd RobotWrapper::joint_positions() const
{
    return utility::to_eigen_vectord(m_robot->jointValues());
}

void RobotWrapper::set_joint_positions(const Eigen::VectorXd &joint_positions)
{
    m_robot->setJointValues(utility::to_std_vectorf(joint_positions));
}
