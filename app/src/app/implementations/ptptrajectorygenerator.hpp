
#ifndef AIS4104_PORTFOLIO_PTPTRAJECTORYGENERATOR_HPP
#define AIS4104_PORTFOLIO_PTPTRAJECTORYGENERATOR_HPP

#include <simulation/trajectorygenerator.h>

namespace AIS4104 {
    // Simple trajectory generator with instantaneous acceleration jumps, Third Order Polynomial
    // Section 9.2.2.1, page 331-332, MR pre-print 2019
    class MyPTPTrajectoryGenerator : public Simulation::PointToPointTrajectoryGenerator
    {
    public:
        using PointToPointTrajectoryGenerator::PointToPointTrajectoryGenerator;

        void stop() override;
        bool has_reached_endpoint() const override;

        bool plan_trajectory(const Simulation::JointLimits &limits, double velocity_factor) override;

        Eigen::VectorXd joint_positions(std::chrono::nanoseconds delta_t) override;

    private:
        std::atomic<bool> m_stopped;
        std::chrono::nanoseconds m_total;
        std::chrono::nanoseconds m_current;
    };

}


#endif //AIS4104_PORTFOLIO_PTPTRAJECTORYGENERATOR_HPP