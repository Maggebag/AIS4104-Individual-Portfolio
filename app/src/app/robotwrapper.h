#ifndef AIS4104_ROBOTWRAPPER_H
#define AIS4104_ROBOTWRAPPER_H

#include <simulation/robot.h>
#include <simulation/kinematicssolver.h>

#include <threepp/objects/Robot.hpp>

namespace AIS4104 {

class RobotWrapper : public Simulation::Robot
{
public:
    RobotWrapper(std::shared_ptr<threepp::Robot> robot, std::shared_ptr<Simulation::KinematicsSolver> solver);

    threepp::Robot& threepp_robot();
    std::shared_ptr<threepp::Robot> threepp_robot_ptr();
    const threepp::Robot& threepp_robot() const;

    uint8_t joint_count() const override;

    Eigen::Matrix4d tool_transform() const override;
    void set_tool_transform(Eigen::Matrix4d transform) override;

    Eigen::Matrix4d current_pose() const override;
    Eigen::Vector3d current_position() const override;
    Eigen::Vector3d current_orientation_zyx() const override;

    Eigen::Matrix4d current_flange_pose() const override;
    Eigen::Vector3d current_flange_position() const override;
    Eigen::Vector3d current_flange_orientation_zyx() const override;

    Eigen::VectorXd ik_solve_pose(const Eigen::Matrix4d &eef_pose, const Eigen::VectorXd &j0) const override;
    Eigen::VectorXd ik_solve_flange_pose(const Eigen::Matrix4d &flange_pose, const Eigen::VectorXd &j0) const override;

    const Simulation::JointLimits& joint_limits() const override;

    Eigen::VectorXd joint_positions() const override;

    void set_joint_positions(const Eigen::VectorXd &joint_positions) override;

private:
    Eigen::Matrix4d m_tool_transform;
    std::shared_ptr<threepp::Robot> m_robot;
    std::shared_ptr<Simulation::KinematicsSolver> m_solver;

};

}

#endif
