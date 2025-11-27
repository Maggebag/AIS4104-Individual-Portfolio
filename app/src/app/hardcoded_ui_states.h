#ifndef AIS4104_HARDCODED_UI_STATES_H
#define AIS4104_HARDCODED_UI_STATES_H

#include <visualization/imgui/windows/toolwindow.h>
#include <visualization/imgui/windows/robotcontrolwindows.h>
#include <visualization/imgui/windows/trajectoryloggerwindow.h>
#include <visualization/imgui/windows/worldobjectconfigurationwindow.h>

namespace AIS4104 {
inline Visualization::TrajectoryLoggerWindow::State default_logger_state()
{
    return Visualization::TrajectoryLoggerWindow::State{
        //Inactive logger by default.
        false,
        //Folder to output logs.
#ifdef _WINDOWS
        "C:/dev/"
#else
        std::string(std::getenv("HOME")) + "/AIS4104"
#endif
    };
}

inline Eigen::VectorXd initial_robot_configuration()
{
    Eigen::VectorXd config(6);
    config << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    return config * utility::deg_to_rad;
}

//Completed-TASK: Default tool configuration:
// a) Set the default value of the orientation of the 3D model (adapter + humerus) to render correctly on the flange.
// b) Set the default value of the frame transformation from the flange to the humeral rotation center.
// (The humerus must face the glenoid for the identity orientation.)

inline Visualization::ToolConfigurationWindow::State default_tool_state()
{
    return Visualization::ToolConfigurationWindow::State{
        //Activate the default tool 3D model (change to show humerus when app opens).
        false,
        //The default tool 3D model to load
        "models/adapter-humerus.stl",
        //The desired view of the tool config window
        Visualization::ToolConfigurationWindow::ToolView::KINEMATICS_TFX,
        //Scale of the tool 3D model used in rendering
        Eigen::Vector3f{1.f, 1.f, 1.f},
        //Position of the tool 3D model used in rendering
        Eigen::Vector3f{0.f, 0.f, 0.f},
        //Euler ZYX orientation (IN DEGREES) of the tool 3D model used in rendering
        Eigen::Vector3f{0.f, -90.f, 0.f},
        //Tool center point offset used for robot kinematics.
        Eigen::Vector3f{0.0425f, -0.025f, 0.1525f},
        // Eigen::Vector3f{0.f, 0.f, 0.f},
        //Euler ZYX orientation (IN DEGREES) of the tool center point used for robot kinematics.
        Eigen::Vector3f{-90.f, -15.f, 0.f}
    };
}

//Keep this as-is
inline Visualization::WorldObjectConfigurationWindow::State default_world_object_state()
{
    return Visualization::WorldObjectConfigurationWindow::State{
        //Activate the default world object 3D model (change to show rig when app opens).
        false,
        //The default world object 3D model to load
        "models/humerus-rig.stl",
        Visualization::WorldObjectConfigurationWindow::WorldView::TRANSFORM,
        //Scale of the world object 3D model used in rendering
        Eigen::Vector3f{0.83f, 0.83f, 0.83f},
        //Position of the world object 3D model used in rendering
        Eigen::Vector3f{1.015f, 0.f, -0.18f},
        //Euler ZYX orientation (IN DEGREES) of the world object 3D model used in rendering
        Eigen::Vector3f{0.f, 90.f, -90.f}
    };
}
}
#endif
