#include "app/implementations/screwskinematicssolver.h"

#include <utility/math.h>

using namespace AIS4104;

ScrewsKinematicsSolver::ScrewsKinematicsSolver(Eigen::Matrix4d m, std::vector<Eigen::VectorXd> screws, Simulation::JointLimits limits)
    : ScrewsKinematicsSolver(std::move(m), std::move(screws), 4.e-3, 4.e-3, std::move(limits))
{
}

ScrewsKinematicsSolver::ScrewsKinematicsSolver(Eigen::Matrix4d m, std::vector<Eigen::VectorXd> space_screws, double v_e, double w_e, Simulation::JointLimits limits)
    : KinematicsSolver(std::move(limits))
    , m_ve(v_e)
    , m_we(w_e)
    , m_m(std::move(m))
    , m_screws(std::move(space_screws))
{
}

void ScrewsKinematicsSolver::set_epsilons(double v_e, double w_e)
{
    m_ve = v_e;
    m_we = w_e;
}

uint32_t ScrewsKinematicsSolver::joint_count() const
{
    return m_screws.size();
}

//TASK: Implement fk_solve using screws.
Eigen::Matrix4d ScrewsKinematicsSolver::fk_solve(const Eigen::VectorXd &joint_positions)
{
    return Eigen::Matrix4d::Identity();
}

Eigen::VectorXd ScrewsKinematicsSolver::ik_solve(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd &j0)
{
    return ik_solve(t_sd, j0, [&](const std::vector<Eigen::VectorXd> &) { return 0u; });
}

//TASK: Implement ik_solve using screws.
Eigen::VectorXd ScrewsKinematicsSolver::ik_solve(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd &j0, const std::function<uint32_t(const std::vector<Eigen::VectorXd> &)> &solution_selector)
{
    return j0;
}

std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ScrewsKinematicsSolver::space_chain()
{
    return {m_m, m_screws};
}

//TASK: Implement body_chain(). You can obtain the variables to transform to body frame from space_chain().
std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ScrewsKinematicsSolver::body_chain()
{
    // auto [m, screws] = space_chain();
    return space_chain();
}

//TASK: Implement space_jacobian() using space_chain()
Eigen::MatrixXd ScrewsKinematicsSolver::space_jacobian(const Eigen::VectorXd &current_joint_positions)
{
    return Eigen::MatrixXd::Identity(current_joint_positions.size(), current_joint_positions.size());
}

//TASK: Implement body_jacobian() using body_chain()
Eigen::MatrixXd ScrewsKinematicsSolver::body_jacobian(const Eigen::VectorXd &current_joint_positions)
{
    return Eigen::MatrixXd::Identity(current_joint_positions.size(), current_joint_positions.size());
}
