#ifndef AIS4104_SIMULATION_ROBOTCONTROLINTERFACE_H
#define AIS4104_SIMULATION_ROBOTCONTROLINTERFACE_H

#include <simulation/robot.h>

#include <Eigen/Dense>

#include <vector>
#include <filesystem>

namespace AIS4104::Simulation {

struct LoggerParameters
{
    bool active = false;
    std::filesystem::path directory = "C:/dev";
};

class RobotControlInterface
{
public:
    RobotControlInterface() = default;
    virtual ~RobotControlInterface() = default;

    virtual LoggerParameters logger_parameters() const = 0;
    virtual void set_logger_parameters(const LoggerParameters &params) = 0;

    virtual uint8_t joint_count() const = 0;

    virtual Eigen::Matrix4d eef_pose() const = 0;
    virtual Eigen::Vector3d eef_position() const = 0;
    virtual Eigen::Vector3f eef_positionf() const = 0;
    virtual Eigen::Vector3d eef_orientation_zyx() const = 0;
    virtual Eigen::Vector3f eef_orientation_zyxf() const = 0;

    virtual Eigen::Matrix4d flange_pose() const = 0;
    virtual Eigen::Vector3d flange_position() const = 0;
    virtual Eigen::Vector3f flange_positionf() const = 0;
    virtual Eigen::Vector3d flange_orientation_zyx() const = 0;
    virtual Eigen::Vector3f flange_orientation_zyxf() const = 0;

    virtual Eigen::VectorXd joint_positions() const = 0;
    virtual Eigen::VectorXf joint_positionsf() const = 0;

    virtual const JointLimits& joint_limits() const = 0;

    virtual double velocity_factor() const = 0;
    virtual void velocity_factor_changed(double factor) = 0; //factor between 0.0 and 1.0 (0-100%)

    virtual uint32_t max_planning_duration() const = 0;
    virtual void max_planning_duration_changed(uint32_t milliseconds) = 0;

    virtual uint32_t max_planning_iterations() const = 0;
    virtual void max_planning_iterations_changed(uint32_t iterations) = 0;

    virtual void on_preview_tool_frame_jog(const Eigen::Matrix4d &start_pose, const Eigen::Vector3d &position, const Eigen::Vector3d &euler_zyx) = 0;
    virtual void on_preview_task_space_pose(const Eigen::Matrix4d &pose) = 0;
    virtual void on_preview_task_space_screw(const Eigen::Matrix4d &start_pose, const Eigen::Vector3d &w, const Eigen::Vector3d &q, double theta, double h) = 0;
    virtual void on_preview_joint_configuration(const Eigen::VectorXd &joint_position) = 0;

    virtual void on_task_space_ptp(const Eigen::Vector3d &pos, const Eigen::Vector3d &e_zyx) = 0;
    virtual void on_task_space_lin(const Eigen::Vector3d &pos, const Eigen::Vector3d &e_zyx) = 0;
    virtual void on_task_space_screw(const Eigen::Vector3d &w, const Eigen::Vector3d &q, double theta, double h) = 0;

    virtual void on_task_space_trajectory(const std::vector<Eigen::Matrix4d> &poses) = 0;
    virtual void on_joint_space_trajectory(const std::vector<Eigen::VectorXd> &joint_positions) = 0;

};

}

#endif
