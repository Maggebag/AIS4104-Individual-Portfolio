#ifndef AIS4104_VISUALIZATION_ROBOTCONTROLWINDOWS_H
#define AIS4104_VISUALIZATION_ROBOTCONTROLWINDOWS_H

#include "visualization/imgui/imguiwindowcontext.h"

#include <simulation/robot.h>
#include <simulation/robotcontrolinterface.h>

#include <utility/enum.h>

namespace AIS4104::Visualization {
enum class ControlMode
{
    PREVIEW, SIMULATION
};

enum class ControlSpace
{
    JOINT, TASK, TOOL
};

enum class Trajectory
{
    PTP, LIN, SCREW
};

class RobotControlWindows : public ImguiWindow
{
    enum class ControlWidgets
    {
        SLIDERS,
        FIELDS
    };

public:
    struct State
    {
        int duration;
        int iterations;
        float velocity;
        Trajectory trajectory;
        ControlMode control_mode;
        ControlSpace control_space;
    };

    RobotControlWindows(std::shared_ptr<Simulation::RobotControlInterface> control_interface);
    RobotControlWindows(std::shared_ptr<Simulation::RobotControlInterface> control_interface, State state);

    State state() const;

    void render() override;

private:
    float m_theta;
    float m_pitch;
    float m_velocity;
    int m_duration;
    int m_iterations;
    Eigen::Vector3f m_q;
    Eigen::Vector3f m_w;
    Eigen::Vector3f m_eef_pos;
    Eigen::Vector3f m_eef_zyx;
    Eigen::Vector3f m_jog_pos_displacement;
    Eigen::Vector3f m_jog_zyx_displacement;
    Eigen::VectorXf m_joint_pos;
    std::vector<std::string> m_joint_labels;
    utility::Enum<Trajectory, 3> m_trajectory;
    utility::Enum<ControlMode, 2> m_control_mode;
    utility::Enum<ControlSpace, 3> m_control_space;
    utility::Enum<ControlWidgets, 2> m_control_widgets;
    std::shared_ptr<Simulation::RobotControlInterface> m_control_interface;

    void render_control_options();
    void render_control_parameters();
    void render_control_inputs();
    void render_joint_space();
    void render_task_space_preview();
    void render_tool_frame_preview();
    void render_task_space_lin_p2p();
    void render_task_space_screw();

    void reset_to_current_pose();

    Eigen::Matrix4d to_matrix(const Eigen::Vector3f &position, const Eigen::Vector3f &euler_zyx);
};
}

#endif
