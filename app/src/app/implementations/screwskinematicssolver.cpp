#include "app/implementations/screwskinematicssolver.h"

#include <iostream>
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

//Completed-TASK: Implement fk_solve using screws.
// Equation 4.14, page 140, MR pre-print 2019
Eigen::Matrix4d ScrewsKinematicsSolver::fk_solve(const Eigen::VectorXd &joint_positions)
{
    Eigen::Matrix4d e_values = Eigen::Matrix4d::Identity();

    if (joint_positions.size() != joint_count()) {
        throw std::invalid_argument("fk_solve: joint positions size does not match joint count");
    }

    for (int i = 0; i < joint_count(); i++) {
        const Eigen::Matrix4d e_i = utility::matrix_exponential(m_screws[i], joint_positions(i));
        e_values *= e_i;
    }

    return e_values * m_m;
}

Eigen::VectorXd ScrewsKinematicsSolver::ik_solve(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd &j0)
{
    return ik_solve(t_sd, j0, [&](const std::vector<Eigen::VectorXd> &) { return 0u; });
}

//Completed-ish-TASK: Implement ik_solve using screws.
// During testing this solver struggles with dynamic changes such as working in the task space with sliders. For single calculations it works
// Numerical IK in body frame (Section 6.2.2, MR preprint 2019) Mainly page 228-229
Eigen::VectorXd ScrewsKinematicsSolver::ik_solve(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd &j0, const std::function<uint32_t(const std::vector<Eigen::VectorXd> &)> &solution_selector)
{
    constexpr size_t max_iter = 1000;
    const double gamma  = 0.1;
    const double lambda = 1e-4;
    size_t iter = 0;

    // IK-solution Canditates
    std::vector<Eigen::VectorXd> candidates;
    {
        // Use current joint positions as initial guess
        Eigen::VectorXd angle_guess = j0;

        for (; iter < max_iter; ++iter) {
            // Compute current end-effector position
            Eigen::Matrix4d T_ee = fk_solve(angle_guess);
            Eigen::Matrix4d T_err = (T_ee.inverse() * t_sd).eval();

            auto [S, theta] = utility::matrix_logarithm(T_err);

            Eigen::VectorXd Vb = S * theta;
            Eigen::Vector3d wb = Vb.head<3>();
            Eigen::Vector3d vb = Vb.tail<3>();

            if (wb.norm() < m_we && vb.norm() < m_ve) {
                candidates.push_back(angle_guess);
                std::cout << iter << std::endl;
                break;
            }

            Eigen::MatrixXd Jb = body_jacobian(angle_guess);

            Eigen::MatrixXd I6 = Eigen::MatrixXd::Identity(6,6);
            Eigen::MatrixXd Jb_pinv =
                Jb.transpose() * (Jb * Jb.transpose() + lambda * lambda * I6).inverse();

            angle_guess += gamma * (Jb_pinv * Vb);
            }
        if (candidates.empty()){candidates.push_back(angle_guess);}
        }

    // Strictly speaking we always currently only compute one solution. It would be beneficial for a solution selector
    // if multiple solutions were computed. Right now it serves no purpose for the screw solver. And since we will not use it further
    // no further solutions will be computed except one.

    uint32_t idx = solution_selector(candidates);
    if (idx >= candidates.size()) {
        throw std::out_of_range("ik_solve: solution_selector returned invalid index");
    }
    std::cout << candidates[idx] << std::endl;
    return candidates[idx];
}


std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ScrewsKinematicsSolver::space_chain()
{
    return {m_m, m_screws};
}

//Completed-TASK: Implement body_chain(). You can obtain the variables to transform to body frame from space_chain().
// Equation 4.16, page 147, MR pre-print 2019
std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ScrewsKinematicsSolver::body_chain()
{
    auto [m, screws] = space_chain();
    const Eigen::Matrix4d m_inv = m.inverse();
    const Eigen::Matrix<double, 6,6> Adminv = utility::adjoint_matrix(m_inv);

    std::vector<Eigen::VectorXd> b_screws(screws.size(), Eigen::VectorXd::Zero(6));
    for (int i = 0 ; i < screws.size(); i++) {
        b_screws[i] = Adminv * screws[i];
    }

    return {m, b_screws};
}

//Completed-TASK: Implement space_jacobian() using space_chain()
// Equation 5.11, page 178, MR pre-print 2019
Eigen::MatrixXd ScrewsKinematicsSolver::space_jacobian(const Eigen::VectorXd &current_joint_positions)
{
    auto [m, screws] = space_chain();
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, screws.size());
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();

    for (int i = 0; i < screws.size(); i++) {
        if (i == 0) {
            J.col(i) = screws[i];
        } else {
            J.col(i) = utility::adjoint_matrix(T) * screws[i];
        }

        T *= utility::matrix_exponential(screws[i], current_joint_positions(i));
    }
    return J;
}

//Completed-TASK: Implement body_jacobian() using body_chain() !! Test if its correct with order of multiplication later
// Equation 5.18, page 183, MR pre-print 2019
Eigen::MatrixXd ScrewsKinematicsSolver::body_jacobian(const Eigen::VectorXd &current_joint_positions)
{
    auto [m, screws] = body_chain();
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, screws.size());
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();

    J.col(screws.size()-1) = screws[screws.size()-1];

    for (int i = screws.size()-2; i >= 0; i--) {
        T *= utility::matrix_exponential(screws[i+1], -current_joint_positions(i+1));

        J.col(i) = utility::adjoint_matrix(T) * screws[i];
    }

    return J;
}
