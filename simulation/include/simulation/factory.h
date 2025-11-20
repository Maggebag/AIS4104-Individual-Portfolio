#ifndef AIS4104_SIMULATION_FACTORY_H
#define AIS4104_SIMULATION_FACTORY_H

#include "simulation/motionplanner.h"
#include "simulation/robotcontrolinterface.h"

namespace AIS4104::Simulation::Factory {

std::shared_ptr<RobotControlInterface> create_control_interface(std::shared_ptr<Robot> robot, std::shared_ptr<KinematicsSolver> solver, std::shared_ptr<MotionPlanner> motion_planner);

}

#endif
