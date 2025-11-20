#include "simulation/kinematicssolver.h"

#include <utility/vectors.h>

using namespace AIS4104::Simulation;

const JointLimits& KinematicsSolver::joint_limits() const
{
    return m_joint_limits;
}

const Parameters& KinematicsSolver::current_solver_parameters() const
{
    return m_solver_params;
}

void KinematicsSolver::set_solver_parameters(Parameters solver_params)
{
    m_solver_params = std::move(solver_params);
}

Eigen::Matrix4d KinematicsSolver::fk_solve(const std::vector<double> &joint_positions)
{
    return fk_solve(utility::to_eigen_vectord(joint_positions));
}

Eigen::VectorXd KinematicsSolver::ik_solve(const Eigen::Matrix4d &desired_tool_pose, const std::vector<double> &j0)
{
    return ik_solve(desired_tool_pose, utility::to_eigen_vectord(j0));
}

Eigen::VectorXd KinematicsSolver::ik_solve(const Eigen::Matrix4d &desired_tool_pose, const std::vector<double> &j0, const std::function<uint32_t(const std::vector<Eigen::VectorXd> &)> &solution_selector)
{
    return ik_solve(desired_tool_pose, utility::to_eigen_vectord(j0), solution_selector);
}
