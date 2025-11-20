#ifndef AIS4104_SIMULATION_KINEMATICSSOLVER_H
#define AIS4104_SIMULATION_KINEMATICSSOLVER_H

#include <Eigen/Dense>

#include <chrono>

namespace AIS4104::Simulation {
struct Parameters
{
    uint32_t flags = 0u;
    uint32_t max_iterations = 1500u;
    uint32_t max_duration_msec = 5u;
    double task_space_epsilon = 1.e-5;
    double joint_space_epsilon = 1.e-15;
};

struct JointLimits
{
    Eigen::VectorXd velocity;
    Eigen::VectorXd acceleration;
    Eigen::VectorXd lower_position;
    Eigen::VectorXd upper_position;
};

class KinematicsSolver
{
public:
    KinematicsSolver(JointLimits limits)
        : m_joint_limits(limits)
    {
    }

    KinematicsSolver(JointLimits limits, Parameters parameters)
        : m_solver_params(parameters)
        , m_joint_limits(limits)
    {
    }

    virtual ~KinematicsSolver() = default;

    virtual uint32_t joint_count() const = 0;

    const JointLimits& joint_limits() const;

    const Parameters& current_solver_parameters() const;

    void set_solver_parameters(Parameters solver_params);

    Eigen::Matrix4d fk_solve(const std::vector<double> &joint_positions);
    virtual Eigen::Matrix4d fk_solve(const Eigen::VectorXd &joint_positions) = 0;

    Eigen::VectorXd ik_solve(const Eigen::Matrix4d &desired_tool_pose, const std::vector<double> &j0);
    virtual Eigen::VectorXd ik_solve(const Eigen::Matrix4d &desired_tool_pose, const Eigen::VectorXd &j0) = 0;

    Eigen::VectorXd ik_solve(const Eigen::Matrix4d &desired_tool_pose, const std::vector<double> &j0, const std::function<uint32_t(const std::vector<Eigen::VectorXd> &)> &solution_idx_selector);
    virtual Eigen::VectorXd ik_solve(const Eigen::Matrix4d &desired_tool_pose, const Eigen::VectorXd &j0, const std::function<uint32_t(const std::vector<Eigen::VectorXd> &)> &solution_idx_selector) = 0;

private:
    Parameters m_solver_params;
    JointLimits m_joint_limits;
};
}

#endif
