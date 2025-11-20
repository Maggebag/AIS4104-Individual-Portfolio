#ifndef AIS4104_SIMULATION_TRAJECTORYGENERATOR_H
#define AIS4104_SIMULATION_TRAJECTORYGENERATOR_H

#include "simulation/kinematicssolver.h"

#include <Eigen/Dense>

#include <chrono>

namespace AIS4104::Simulation {

class TrajectoryGenerator
{
public:
    TrajectoryGenerator() = default;
    virtual ~TrajectoryGenerator() = default;

    virtual void stop() = 0;
    virtual bool has_reached_endpoint() const = 0;

    virtual bool plan_trajectory(const JointLimits &limits, double velocity_factor) = 0;

    virtual Eigen::VectorXd joint_positions(std::chrono::nanoseconds delta_t) = 0;
};

class PointToPointTrajectoryGenerator : public TrajectoryGenerator
{
public:
    PointToPointTrajectoryGenerator(const Eigen::VectorXd &w0, const Eigen::VectorXd &w1)
        : m_w0(w0)
        , m_w1(w1)
    {
    }

    ~PointToPointTrajectoryGenerator() override = default;

protected:
    Eigen::VectorXd m_w0;
    Eigen::VectorXd m_w1;
};

class MultipointTrajectoryGenerator : public TrajectoryGenerator
{
public:
    MultipointTrajectoryGenerator(const std::vector<Eigen::VectorXd> &waypoints)
        : m_waypoints(waypoints)
    {
    }

    ~MultipointTrajectoryGenerator() override = default;

protected:
    std::vector<Eigen::VectorXd> m_waypoints;
};

}

#endif
