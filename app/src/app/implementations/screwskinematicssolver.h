#ifndef AIS4104_SCREWKINEMATICSSOLVER_H
#define AIS4104_SCREWKINEMATICSSOLVER_H

#include <simulation/kinematicssolver.h>

#include <chrono>

namespace AIS4104 {

class ScrewsKinematicsSolver : public Simulation::KinematicsSolver
{
public:
    ScrewsKinematicsSolver(Eigen::Matrix4d m, std::vector<Eigen::VectorXd> space_screws, Simulation::JointLimits limits);
    ScrewsKinematicsSolver(Eigen::Matrix4d m, std::vector<Eigen::VectorXd> space_screws, double v_e, double w_e, Simulation::JointLimits limits);

    void set_epsilons(double v_e, double w_e);

    uint32_t joint_count() const override;

    std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> space_chain();
    std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> body_chain();

    Eigen::MatrixXd space_jacobian(const Eigen::VectorXd &current_joint_positions);
    Eigen::MatrixXd body_jacobian(const Eigen::VectorXd &current_joint_positions);

private:
    double m_ve;
    double m_we;
    Eigen::Matrix4d m_m;
    std::vector<Eigen::VectorXd> m_screws;

    Eigen::Matrix4d fk_solve(const Eigen::VectorXd &joint_positions) override;
    Eigen::VectorXd ik_solve(const Eigen::Matrix4d &desired_flange_pose, const Eigen::VectorXd &j0) override;
    Eigen::VectorXd ik_solve(const Eigen::Matrix4d &desired_flange_pose, const Eigen::VectorXd &j0, const std::function<uint32_t(const std::vector<Eigen::VectorXd> &)> &solution_selector) override;
};

}

#endif
