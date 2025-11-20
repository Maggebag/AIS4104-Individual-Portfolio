#ifndef AIS4104_SIMULATION_ROBOTCONTROLDELEGATOR_H
#define AIS4104_SIMULATION_ROBOTCONTROLDELEGATOR_H

#include "simulation/robot.h"
#include "simulation/motionplanner.h"
#include "simulation/kinematicssolver.h"
#include "simulation/trajectoryexecutor.h"
#include "simulation/robotcontrolinterface.h"

namespace AIS4104::Simulation {

class RobotControlDelegator : public RobotControlInterface
{
public:
    RobotControlDelegator(std::shared_ptr<Robot> robot, std::shared_ptr<KinematicsSolver> kinematics_solver, std::shared_ptr<MotionPlanner> motion_planner);

    LoggerParameters logger_parameters() const override;
    void set_logger_parameters(const LoggerParameters &params) override;

    uint8_t joint_count() const override;

    Eigen::Matrix4d eef_pose() const override;
    Eigen::Vector3d eef_position() const override;
    Eigen::Vector3f eef_positionf() const override;
    Eigen::Vector3d eef_orientation_zyx() const override;
    Eigen::Vector3f eef_orientation_zyxf() const override;

    Eigen::Matrix4d flange_pose() const override;
    Eigen::Vector3d flange_position() const override;
    Eigen::Vector3f flange_positionf() const override;
    Eigen::Vector3d flange_orientation_zyx() const override;
    Eigen::Vector3f flange_orientation_zyxf() const override;

    Eigen::VectorXd joint_positions() const override;
    Eigen::VectorXf joint_positionsf() const override;

    const JointLimits& joint_limits() const override;

    double velocity_factor() const override;
    void velocity_factor_changed(double factor) override;

    uint32_t max_planning_duration() const override;
    void max_planning_duration_changed(uint32_t milliseconds) override;

    uint32_t max_planning_iterations() const override;
    void max_planning_iterations_changed(uint32_t iterations) override;

    void on_preview_tool_frame_jog(const Eigen::Matrix4d &start_pose, const Eigen::Vector3d &displacement, const Eigen::Vector3d &rotation) override;
    void on_preview_task_space_pose(const Eigen::Matrix4d &pose) override;
    void on_preview_task_space_screw(const Eigen::Matrix4d &start_pose, const Eigen::Vector3d &w, const Eigen::Vector3d &q, double theta, double h) override;
    void on_preview_joint_configuration(const Eigen::VectorXd &joint_position) override;

    void on_task_space_ptp(const Eigen::Vector3d &pos, const Eigen::Vector3d &e_zyx) override;
    void on_task_space_lin(const Eigen::Vector3d &pos, const Eigen::Vector3d &e_zyx) override;
    void on_task_space_screw(const Eigen::Vector3d &w, const Eigen::Vector3d &q, double theta, double h) override;

    void on_task_space_trajectory(const std::vector<Eigen::Matrix4d> &poses) override;
    void on_joint_space_trajectory(const std::vector<Eigen::VectorXd> &joint_positions) override;

private:
    std::shared_ptr<Robot> m_robot;
    std::shared_ptr<MotionPlanner> m_motion_planner;
    std::shared_ptr<KinematicsSolver> m_kinematics_solver;
    std::shared_ptr<TrajectoryExecutor> m_trajectory_executor;

    void run_trajectory(const std::shared_ptr<TrajectoryGenerator> &generator);
};

}

#endif
