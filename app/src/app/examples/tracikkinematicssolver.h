#ifndef AIS4104_APP_TRACIKKINEMATICSSOLVER_H
#define AIS4104_APP_TRACIKKINEMATICSSOLVER_H

#include <simulation/kinematicssolver.h>

#include <Eigen/Dense>

#include <trac_ik/trac_ik.hpp>

#include <chrono>

namespace AIS4104 {

class TracIkKinematicsSolver : public Simulation::KinematicsSolver
{
public:
    TracIkKinematicsSolver(std::shared_ptr<KDL::Chain> kinematic_chain, Simulation::JointLimits limits);

    uint32_t joint_count() const override;

    Eigen::Matrix4d fk_solve(const Eigen::VectorXd &joint_positions) override;

    Eigen::VectorXd ik_solve(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd &j0) override;
    Eigen::VectorXd ik_solve(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd &j0, const std::function<uint32_t(const std::vector<Eigen::VectorXd> &)> &solution_selector) override;

private:
    KDL::JntArray m_lower_joint_limits;
    KDL::JntArray m_upper_joint_limits;
    Simulation::Parameters m_solver_params;
    std::shared_ptr<KDL::Chain> m_kinematic_chain;

};

}

#endif
