
#include "ptptrajectorygenerator.hpp"

using namespace AIS4104;

void MyPTPTrajectoryGenerator::stop()
{
    m_stopped = true;
}

bool MyPTPTrajectoryGenerator::has_reached_endpoint() const
{
    return m_current == m_total || m_stopped;
}

bool MyPTPTrajectoryGenerator::plan_trajectory(const Simulation::JointLimits &limits, double velocity_factor)
{
}

Eigen::VectorXd MyPTPTrajectoryGenerator::joint_positions(std::chrono::nanoseconds delta_t)
{
}
