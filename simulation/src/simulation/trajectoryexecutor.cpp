#include "simulation/trajectoryexecutor.h"

#include <utility/csv.h>

#include "spdlog/spdlog.h"

using namespace AIS4104::Simulation;

TrajectoryExecutor::TrajectoryExecutor(std::shared_ptr<Robot> robot)
    : m_do_work{true}
    , m_executing{false}
    , m_robot(std::move(robot))
{
    m_execution_worker = std::thread([this]()
        {
            while(m_do_work)
            {
                std::unique_lock l(m_trajectory_mutex);
                m_execution_condition.wait(l, [this]() { return is_executing_trajectory(); });
                m_executing = true;
                unsafe_execute_trajectory();
            }
        }
    );
}

TrajectoryExecutor::~TrajectoryExecutor()
{
    m_do_work = false;
    if(m_executing)
        stop();
    if(m_execution_worker.joinable())
        m_execution_worker.join();
}

void TrajectoryExecutor::set_logger_parameters(LoggerParameters logger_parameters)
{
    std::lock_guard<std::mutex> log_lock(m_log_mutex);
    m_log_params = std::move(logger_parameters);
}

LoggerParameters TrajectoryExecutor::logger_parameters() const
{
    std::lock_guard<std::mutex> log_lock(m_log_mutex);
    return m_log_params;
}

void TrajectoryExecutor::stop()
{
    m_executing = false;
}

bool TrajectoryExecutor::is_executing_trajectory() const
{
    return m_executing;
}

bool TrajectoryExecutor::execute_trajectory(std::shared_ptr<TrajectoryGenerator> trajectory)
{
    if(trajectory == nullptr)
        return false;
    {
        std::unique_lock<std::mutex> traj_l(m_trajectory_mutex, std::try_to_lock);
        if(!traj_l.owns_lock())
            return false;
        m_trajectory = std::move(trajectory);
        m_executing = true;
    }
    m_execution_condition.notify_one();
    return true;
}

void TrajectoryExecutor::unsafe_execute_trajectory()
{
    if(!m_trajectory->plan_trajectory(m_robot->joint_limits(), m_robot->velocity_factor()))
    {
        m_executing = false;
        spdlog::error("[TrajectoryGenerator::plan_trajectory] Unknown error occurred: unable to plan trajectory.");
    }
    LoggerParameters log_params;
    {
        std::lock_guard<std::mutex> log_lock(m_log_mutex);
        log_params = m_log_params;
    }
    std::chrono::time_point<std::chrono::steady_clock> previous_time = std::chrono::steady_clock::now();
    std::vector<uint64_t> timestamps;
    std::vector<Eigen::VectorXd> trajectory;
    while(m_executing && !m_trajectory->has_reached_endpoint())
    {
        auto now = std::chrono::steady_clock::now();
        auto delta_t = std::chrono::duration_cast<std::chrono::nanoseconds>(now - previous_time);
        Eigen::VectorXd jpos = m_trajectory->joint_positions(delta_t);
        m_robot->set_joint_positions(jpos);
        previous_time = now;
        if(log_params.active)
        {
            trajectory.push_back(jpos);
            timestamps.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count());
        }
    }
    if(log_params.active)
    {
        utility::std_vector_to_csv_file(timestamps, log_params.directory / utility::date_time_stamped_filename("timestamp.csv"));
        utility::eigen_vector_to_csv_file(trajectory, log_params.directory / utility::date_time_stamped_filename("trajectory.csv"));
    }
    m_executing = false;
}
