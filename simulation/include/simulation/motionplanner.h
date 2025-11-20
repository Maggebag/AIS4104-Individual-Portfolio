#ifndef AIS4104_SIMULATION_MOTIONPLANNER_H
#define AIS4104_SIMULATION_MOTIONPLANNER_H

#include "simulation/trajectorygenerator.h"

#include <memory>

namespace AIS4104::Simulation {

class MotionPlanner
{
public:
    MotionPlanner() = default;
    virtual ~MotionPlanner() = default;

    virtual Eigen::VectorXd task_space_pose(const Eigen::Matrix4d &pose) = 0;
    virtual Eigen::VectorXd task_space_screw(const Eigen::Matrix4d &world_start_pose, const Eigen::Vector3d &w, const Eigen::Vector3d &q, double theta, double h) = 0;
    virtual Eigen::VectorXd tool_frame_displace(const Eigen::Matrix4d &world_start_pose, const Eigen::Vector3d &tf_offset, const Eigen::Vector3d &tf_zyx) = 0;

    virtual std::shared_ptr<TrajectoryGenerator> task_space_ptp_trajectory(const Eigen::Vector3d &target_pos, const Eigen::Vector3d &target_e_zyx) = 0;
    virtual std::shared_ptr<TrajectoryGenerator> task_space_lin_trajectory(const Eigen::Vector3d &target_pos, const Eigen::Vector3d &target_e_zyx) = 0;
    virtual std::shared_ptr<TrajectoryGenerator> task_space_screw_trajectory(const Eigen::Vector3d &w, const Eigen::Vector3d &q, double theta, double h) = 0;

    virtual std::shared_ptr<TrajectoryGenerator> task_space_trajectory(const std::vector<Eigen::Matrix4d> &waypoints) = 0;
    virtual std::shared_ptr<TrajectoryGenerator> joint_space_trajectory(const std::vector<Eigen::VectorXd> &waypoints) = 0;
};

}

#endif
