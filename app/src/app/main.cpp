#include "app/robotwrapper.h"

#include "app/urdf_loaders.h"
#include "app/hardcoded_loaders.h"
#include "app/hardcoded_ui_states.h"

#include "app/implementations/motionplannerimpl.h"

#include <simulation/factory.h>

#include <visualization/robotscene.h>

#include <visualization/imgui/windows/toolwindow.h>
#include <visualization/imgui/windows/posedisplaywindow.h>
#include <visualization/imgui/windows/robotcontrolwindows.h>
#include <visualization/imgui/windows/trajectoryloggerwindow.h>
#include <visualization/imgui/windows/worldobjectconfigurationwindow.h>

void render_scene_without_default_tool(std::shared_ptr<AIS4104::RobotWrapper> robot, std::shared_ptr<AIS4104::Simulation::RobotControlInterface> control_iface)
{
    AIS4104::Visualization::RobotScene scene(robot->threepp_robot_ptr(), {
            std::make_shared<AIS4104::Visualization::WorldObjectConfigurationWindow>(scene),
            std::make_shared<AIS4104::Visualization::RobotControlWindows>(control_iface),
            std::make_shared<AIS4104::Visualization::PoseDisplayWindow>(*control_iface),
            std::make_shared<AIS4104::Visualization::ToolConfigurationWindow>(scene, robot),
            std::make_shared<AIS4104::Visualization::TrajectoryLoggerWindow>(*control_iface)
        },
        AIS4104::Visualization::RobotScene::Projection::ORTHOGRAPHIC
    );
    scene.render();
}

void render_scene_with_default_tool(std::shared_ptr<AIS4104::RobotWrapper> robot, std::shared_ptr<AIS4104::Simulation::RobotControlInterface> control_iface)
{
    AIS4104::Visualization::RobotScene scene(robot->threepp_robot_ptr(), {
            std::make_shared<AIS4104::Visualization::WorldObjectConfigurationWindow>(scene, AIS4104::default_world_object_state()),
            std::make_shared<AIS4104::Visualization::RobotControlWindows>(control_iface),
            std::make_shared<AIS4104::Visualization::PoseDisplayWindow>(*control_iface),
            std::make_shared<AIS4104::Visualization::ToolConfigurationWindow>(scene, robot, AIS4104::default_tool_state()),
            std::make_shared<AIS4104::Visualization::TrajectoryLoggerWindow>(*control_iface, AIS4104::default_logger_state())
        },
        AIS4104::Visualization::RobotScene::Projection::ORTHOGRAPHIC
    );
    scene.show_grid(true);
    scene.set_background_color(threepp::Color::white);
    scene.render();
}

int main(int, char **)
{
    std::filesystem::path kr6_path = "urdf/urdf_files/ros-industrial/kuka_kr6_support/urdf/kr6r900sixx.urdf";
    std::filesystem::path ur3e_path = "urdf/urdf_files/matlab/ur_e_description/urdf/universalUR3e.urdf";

    auto solver = AIS4104::hardcoded_kr6r_tracik_solver();
    auto robot = std::make_shared<AIS4104::RobotWrapper>(AIS4104::load_threepp_robot(kr6_path), solver);
    auto control_iface = AIS4104::Simulation::Factory::create_control_interface(robot, solver, std::make_shared<AIS4104::MotionPlannerImpl>(*robot));
    robot->set_joint_positions(AIS4104::initial_robot_configuration());

    // render_scene_without_default_tool(robot, control_iface);
    render_scene_with_default_tool(robot, control_iface);

    return 0;
}
