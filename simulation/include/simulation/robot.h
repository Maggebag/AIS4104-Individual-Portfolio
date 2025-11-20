#ifndef AIS4104_SIMULATION_ROBOT_H
#define AIS4104_SIMULATION_ROBOT_H

#include "simulation/kinematicssolver.h"

#include <Eigen/Dense>

namespace AIS4104::Simulation {

class Robot
{
public:
    Robot()
        : m_velocity_factor(0.3)
    {
    }

    virtual ~Robot() = default;

    double velocity_factor() const
    {
        return m_velocity_factor;
    }

    void set_velocity_factor(double factor)
    {
        m_velocity_factor = factor;
    }

    virtual uint8_t joint_count() const = 0;

    virtual Eigen::Matrix4d tool_transform() const = 0;
    virtual void set_tool_transform(Eigen::Matrix4d transform) = 0;

    virtual Eigen::Matrix4d current_pose() const = 0;
    virtual Eigen::Vector3d current_position() const = 0;
    virtual Eigen::Vector3d current_orientation_zyx() const = 0;

    virtual Eigen::Matrix4d current_flange_pose() const = 0;
    virtual Eigen::Vector3d current_flange_position() const = 0;
    virtual Eigen::Vector3d current_flange_orientation_zyx() const = 0;

    virtual Eigen::VectorXd ik_solve_pose(const Eigen::Matrix4d &tf, const Eigen::VectorXd &j0) const = 0;
    virtual Eigen::VectorXd ik_solve_flange_pose(const Eigen::Matrix4d &tf, const Eigen::VectorXd &j0) const = 0;

    virtual Eigen::VectorXd joint_positions() const = 0;
    virtual void set_joint_positions(const Eigen::VectorXd &joint_positions) = 0;

    virtual const JointLimits& joint_limits() const = 0;

private:
    std::atomic<double> m_velocity_factor;
};

}

#endif
