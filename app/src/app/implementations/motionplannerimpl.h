#ifndef AIS4104_APP_ROBOTCONTROLDELEGATOR_H
#define AIS4104_APP_ROBOTCONTROLDELEGATOR_H

#include <simulation/robot.h>
#include <simulation/motionplanner.h>

namespace AIS4104 {

class MotionPlannerImpl : public Simulation::MotionPlanner
{
public:
    MotionPlannerImpl(const Simulation::Robot &robot);

    Eigen::VectorXd task_space_pose(const Eigen::Matrix4d &pose) override;
    Eigen::VectorXd task_space_screw(const Eigen::Matrix4d &tw_start_pose, const Eigen::Vector3d &w, const Eigen::Vector3d &q, double theta, double h) override;
    Eigen::VectorXd tool_frame_displace(const Eigen::Matrix4d &tw_start_pose, const Eigen::Vector3d &tf_offset, const Eigen::Vector3d &tf_zyx) override;

    std::shared_ptr<Simulation::TrajectoryGenerator> task_space_ptp_trajectory(const Eigen::Vector3d &pos, const Eigen::Vector3d &euler_zyx) override;
    std::shared_ptr<Simulation::TrajectoryGenerator> task_space_lin_trajectory(const Eigen::Vector3d &pos, const Eigen::Vector3d &euler_zyx) override;
    std::shared_ptr<Simulation::TrajectoryGenerator> task_space_screw_trajectory(const Eigen::Vector3d &w, const Eigen::Vector3d &q, double theta, double h) override;

    std::shared_ptr<Simulation::TrajectoryGenerator> task_space_trajectory(const std::vector<Eigen::Matrix4d> &waypoints) override;
    std::shared_ptr<Simulation::TrajectoryGenerator> joint_space_trajectory(const std::vector<Eigen::VectorXd> &waypoints) override;

private:
    const Simulation::Robot &m_robot;
};

}

#endif
