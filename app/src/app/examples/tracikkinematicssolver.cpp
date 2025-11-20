#include "app/examples/tracikkinematicssolver.h"

#include <utility/math.h>
#include <utility/vectors.h>

#include <trac_ik/trac_ik.hpp>

#include "spdlog/spdlog.h"

using namespace AIS4104;

TracIkKinematicsSolver::TracIkKinematicsSolver(std::shared_ptr<KDL::Chain> kinematic_chain, Simulation::JointLimits limits)
: KinematicsSolver(limits)
, m_lower_joint_limits(utility::to_kdl_array(joint_limits().lower_position))
, m_upper_joint_limits(utility::to_kdl_array(joint_limits().upper_position))
, m_kinematic_chain(std::move(kinematic_chain))
{
    m_solver_params.flags = static_cast<uint32_t>(TRAC_IK::Manip2);
}

uint32_t TracIkKinematicsSolver::joint_count() const
{
    return m_kinematic_chain->getNrOfJoints();
}

Eigen::Matrix4d TracIkKinematicsSolver::fk_solve(const Eigen::VectorXd &joint_positions)
{
    KDL::ChainFkSolverPos_recursive solver(*m_kinematic_chain);
    KDL::Frame out;
    solver.JntToCart(utility::to_kdl_array(joint_positions), out);
    return utility::transformation_matrix(utility::to_eigen_matrix3d(out.M), utility::to_eigen_vector3d(out.p));
}

Eigen::VectorXd TracIkKinematicsSolver::ik_solve(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd &j0)
{
    return ik_solve(t_sd, j0, [&](const std::vector<Eigen::VectorXd> &v) { return 0u; });
}

Eigen::VectorXd TracIkKinematicsSolver::ik_solve(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd &j0, const std::function<uint32_t(const std::vector<Eigen::VectorXd> &)> &solution_selector)
{
    TRAC_IK::TRAC_IK solver(*m_kinematic_chain, m_lower_joint_limits, m_upper_joint_limits, m_solver_params.max_duration_msec / 1000.0, m_solver_params.task_space_epsilon, static_cast<TRAC_IK::SolveType>(m_solver_params.flags));
    KDL::JntArray solution;
    std::vector<KDL::JntArray> kdl_solutions;
    auto result = solver.CartToJnt(utility::to_kdl_array(j0), utility::to_kdl_frame(t_sd), solution);
    if(result > 0)
    {
        solver.getSolutions(kdl_solutions);
        std::vector<Eigen::VectorXd> eigen_solutions(kdl_solutions.size());
        for(auto i = 0u; i < kdl_solutions.size(); ++i)
            eigen_solutions[i] = std::move(kdl_solutions[i].data);
        auto jt = eigen_solutions[solution_selector(eigen_solutions)];
        if(jt.hasNaN())
        {
            spdlog::error("IK solution has NaN!");
            return j0;
        }
        return jt;
    }
    return j0;
}
