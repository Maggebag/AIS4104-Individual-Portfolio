#ifndef AIS4104_SIMULATION_TRAJECTORYEXECUTOR_H
#define AIS4104_SIMULATION_TRAJECTORYEXECUTOR_H

#include "simulation/robot.h"
#include "simulation/trajectorygenerator.h"
#include "simulation/robotcontrolinterface.h"

#include <mutex>
#include <atomic>
#include <thread>
#include <condition_variable>

namespace AIS4104::Simulation {

class TrajectoryExecutor
{
public:
    TrajectoryExecutor(std::shared_ptr<Robot> robot);
    ~TrajectoryExecutor();

    void set_logger_parameters(LoggerParameters logger_parameters);
    LoggerParameters logger_parameters() const;

    void stop();

    bool is_executing_trajectory() const;

    bool execute_trajectory(std::shared_ptr<TrajectoryGenerator> trajectory);

private:
    std::mutex m_trajectory_mutex;
    mutable std::mutex m_log_mutex;
    std::thread m_execution_worker;
    LoggerParameters m_log_params;
    std::atomic<bool> m_do_work;
    std::atomic<bool> m_executing;
    std::condition_variable m_execution_condition;
    std::shared_ptr<Robot> m_robot;
    std::shared_ptr<TrajectoryGenerator> m_trajectory;

    void unsafe_execute_trajectory();
};

}

#endif
