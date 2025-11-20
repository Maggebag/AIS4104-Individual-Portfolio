#include "app/examples/ptplinetrajectorygeneratorexample.h"

using namespace AIS4104;

void PTPLineTrajectoryGeneratorExample::stop()
{
    m_stopped = true;
}

bool PTPLineTrajectoryGeneratorExample::has_reached_endpoint() const
{
    return m_current == m_total || m_stopped;
}

bool PTPLineTrajectoryGeneratorExample::plan_trajectory(const Simulation::JointLimits &limits, double velocity_factor)
{
    if(m_w0.isApprox(m_w1))
    {
        m_total = std::chrono::nanoseconds(0);
        m_stopped = true;
    }
    else
    {
        m_stopped = false;
        Eigen::VectorXd velocity_limits = limits.velocity * velocity_factor;
        Eigen::VectorXd diff = (m_w1 - m_w0).cwiseAbs();
        Eigen::VectorXd times = diff.array() / velocity_limits.array();
        uint64_t duration = (times.maxCoeff() * 1000000000.0);
        m_total = std::chrono::nanoseconds(duration);
    }
    return true;
}

Eigen::VectorXd PTPLineTrajectoryGeneratorExample::joint_positions(std::chrono::nanoseconds delta_t)
{
    if(m_stopped)
    {
        double s = std::clamp(static_cast<double>(m_current.count()) / static_cast<double>(m_total.count()), 0.0, 1.0);
        return m_w0 + s * (m_w1 - m_w0);
    }
    m_current = std::clamp(m_current + delta_t, std::chrono::nanoseconds(0), m_total);
    double s = std::clamp(static_cast<double>(m_current.count()) / static_cast<double>(m_total.count()), 0.0, 1.0);
    return m_w0 + s * (m_w1 - m_w0);
}
