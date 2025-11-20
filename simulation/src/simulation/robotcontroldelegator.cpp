#include "simulation/robotcontroldelegator.h"

#include "simulation/trajectorygenerator.h"

#include <utility/vectors.h>

#include "spdlog/spdlog.h"

using namespace AIS4104::Simulation;

RobotControlDelegator::RobotControlDelegator(std::shared_ptr<Robot> robot, std::shared_ptr<KinematicsSolver> kinematics_solver, std::shared_ptr<MotionPlanner> motion_planner)
    : m_robot(std::move(robot))
    , m_motion_planner(std::move(motion_planner))
    , m_kinematics_solver(std::move(kinematics_solver))
    , m_trajectory_executor(std::make_shared<TrajectoryExecutor>(m_robot))
{
}

LoggerParameters RobotControlDelegator::logger_parameters() const
{
    return m_trajectory_executor->logger_parameters();
}

void RobotControlDelegator::set_logger_parameters(const LoggerParameters &params)
{
    m_trajectory_executor->set_logger_parameters(params);
}

uint8_t RobotControlDelegator::joint_count() const
{
    return m_robot->joint_count();
}

Eigen::Matrix4d RobotControlDelegator::eef_pose() const
{
    return m_robot->current_pose();
}

Eigen::Vector3d RobotControlDelegator::eef_position() const
{
    return m_robot->current_position();
}

Eigen::Vector3f RobotControlDelegator::eef_positionf() const
{
    return utility::to_eigen_vectorf(m_robot->current_position());
}

Eigen::Vector3d RobotControlDelegator::eef_orientation_zyx() const
{
    return m_robot->current_orientation_zyx();
}

Eigen::Vector3f RobotControlDelegator::eef_orientation_zyxf() const
{
    return utility::to_eigen_vectorf(m_robot->current_orientation_zyx());
}

Eigen::Matrix4d RobotControlDelegator::flange_pose() const
{
    return m_robot->current_flange_pose();
}

Eigen::Vector3d RobotControlDelegator::flange_position() const
{
    return m_robot->current_flange_position();
}

Eigen::Vector3f RobotControlDelegator::flange_positionf() const
{
    return utility::to_eigen_vectorf(m_robot->current_flange_position());
}

Eigen::Vector3d RobotControlDelegator::flange_orientation_zyx() const
{
    return m_robot->current_flange_orientation_zyx();
}

Eigen::Vector3f RobotControlDelegator::flange_orientation_zyxf() const
{
    return utility::to_eigen_vectorf(m_robot->current_flange_orientation_zyx());
}

Eigen::VectorXd RobotControlDelegator::joint_positions() const
{
    return m_robot->joint_positions();
}

Eigen::VectorXf RobotControlDelegator::joint_positionsf() const
{
    return utility::to_eigen_vectorf(m_robot->joint_positions());
}

const JointLimits& RobotControlDelegator::joint_limits() const
{
    return m_robot->joint_limits();
}

double RobotControlDelegator::velocity_factor() const
{
    return m_robot->velocity_factor();
}

void RobotControlDelegator::velocity_factor_changed(double factor)
{
    m_robot->set_velocity_factor(factor);
}

uint32_t RobotControlDelegator::max_planning_duration() const
{
    return m_kinematics_solver->current_solver_parameters().max_duration_msec;
}

void RobotControlDelegator::max_planning_duration_changed(uint32_t milliseconds)
{
    auto params = m_kinematics_solver->current_solver_parameters();
    params.max_duration_msec = milliseconds;
    m_kinematics_solver->set_solver_parameters(params);
}

uint32_t RobotControlDelegator::max_planning_iterations() const
{
    return m_kinematics_solver->current_solver_parameters().max_iterations;
}

void RobotControlDelegator::max_planning_iterations_changed(uint32_t iterations)
{
    auto params = m_kinematics_solver->current_solver_parameters();
    params.max_duration_msec = iterations;
    m_kinematics_solver->set_solver_parameters(params);
}

void RobotControlDelegator::on_preview_tool_frame_jog(const Eigen::Matrix4d &start_pose, const Eigen::Vector3d &displacement, const Eigen::Vector3d &rotation)
{
    if(m_trajectory_executor->is_executing_trajectory())
        return;
    m_robot->set_joint_positions(m_motion_planner->tool_frame_displace(start_pose, displacement, rotation));
}

void RobotControlDelegator::on_preview_task_space_pose(const Eigen::Matrix4d &pose)
{
    if(m_trajectory_executor->is_executing_trajectory())
        return;
    m_robot->set_joint_positions(m_motion_planner->task_space_pose(pose));
}

void RobotControlDelegator::on_preview_joint_configuration(const Eigen::VectorXd &joint_position)
{
    if(m_trajectory_executor->is_executing_trajectory())
        return;
    m_robot->set_joint_positions(joint_position);
}

void RobotControlDelegator::on_preview_task_space_screw(const Eigen::Matrix4d &start_pose, const Eigen::Vector3d &w, const Eigen::Vector3d &q, double theta, double h)
{
    if(m_trajectory_executor->is_executing_trajectory())
        return;
    auto generator = m_motion_planner->joint_space_trajectory({joint_positions(), m_motion_planner->task_space_screw(start_pose, w, q, theta, h)});
    run_trajectory(generator);
}

void RobotControlDelegator::on_task_space_ptp(const Eigen::Vector3d &pos, const Eigen::Vector3d &e_zyx)
{
    if(m_trajectory_executor->is_executing_trajectory())
        return;
    run_trajectory(m_motion_planner->task_space_ptp_trajectory(pos, e_zyx));
}

void RobotControlDelegator::on_task_space_lin(const Eigen::Vector3d &pos, const Eigen::Vector3d &e_zyx)
{
    if(m_trajectory_executor->is_executing_trajectory())
        return;
    run_trajectory(m_motion_planner->task_space_lin_trajectory(pos, e_zyx));
}

void RobotControlDelegator::on_task_space_screw(const Eigen::Vector3d &w, const Eigen::Vector3d &q, double theta, double h)
{
    run_trajectory(m_motion_planner->task_space_screw_trajectory(w, q, theta, h));
}

void RobotControlDelegator::on_task_space_trajectory(const std::vector<Eigen::Matrix4d> &poses)
{
    run_trajectory(m_motion_planner->task_space_trajectory(poses));
}

void RobotControlDelegator::on_joint_space_trajectory(const std::vector<Eigen::VectorXd> &joint_positions)
{
    run_trajectory(m_motion_planner->joint_space_trajectory(joint_positions));
}

void RobotControlDelegator::run_trajectory(const std::shared_ptr<TrajectoryGenerator> &generator)
{
    if(generator != nullptr)
        m_trajectory_executor->execute_trajectory(generator);
    else
        spdlog::warn("[RobotControlDelegator::run_trajectory] No trajectory executed, nullptr received.");
}
