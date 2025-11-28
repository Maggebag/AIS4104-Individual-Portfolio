
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

// Section 9.2.2.1, page 331-332, MR pre-print 2019
bool MyPTPTrajectoryGenerator::plan_trajectory(const Simulation::JointLimits &limits, double velocity_factor)
{
    if(m_w0.isApprox(m_w1))
    {
        m_total = std::chrono::nanoseconds(0);
        m_current = std::chrono::nanoseconds(0);
        m_stopped = true;
    }
    else
    {
        // Computes the amount of time for the planed trajectory
        m_stopped = false;
        Eigen::VectorXd velocity_limits = limits.velocity * velocity_factor;
        Eigen::VectorXd diff = (m_w1 - m_w0).cwiseAbs();
        Eigen::VectorXd times = 1.5 * diff.array() / velocity_limits.array(); // Added 1.5 times due to polynomial time scaling
        uint64_t duration = (times.maxCoeff() * 1000000000.0);
        m_total = std::chrono::nanoseconds(duration);
    }
    return true;
}

// Section 9.2.2.1, page 331-332, MR pre-print 2019
Eigen::VectorXd MyPTPTrajectoryGenerator::joint_positions(std::chrono::nanoseconds delta_t)
{
    double t = static_cast<double>(m_current.count());
    double T = static_cast<double>(m_total.count());

    if(m_stopped)
    {
        double s = std::clamp(3.0 * t * t / (T * T) - 2.0 * t * t * t / (T * T * T), 0.0, 1.0);
        return m_w0 + s * (m_w1 - m_w0);
    }

    m_current = std::clamp(m_current + delta_t, std::chrono::nanoseconds(0), m_total);
    double s = std::clamp(3.0 * t * t / (T * T) - 2.0 * t * t * t / (T * T * T), 0.0, 1.0);
    return m_w0 + s * (m_w1 - m_w0);
}
