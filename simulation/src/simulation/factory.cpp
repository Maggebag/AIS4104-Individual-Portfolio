#include "simulation/factory.h"

#include "trajectoryexecutor.h"
#include "robotcontroldelegator.h"

namespace AIS4104::Simulation::Factory {

std::shared_ptr<RobotControlInterface> create_control_interface(std::shared_ptr<Robot> robot, std::shared_ptr<KinematicsSolver> solver, std::shared_ptr<MotionPlanner> motion_planner)
{
    return std::make_shared<RobotControlDelegator>(robot, solver, motion_planner);
}

}
