#include "visualization/imgui/windows/robotcontrolwindows.h"

#include "visualization/imgui/widgets.h"

#include <utility/math.h>
#include <utility/vectors.h>

using namespace AIS4104::Visualization;

namespace tpp = threepp;

RobotControlWindows::RobotControlWindows(std::shared_ptr<Simulation::RobotControlInterface> control_interface)
    : m_theta(0.f)
    , m_pitch(0.f)
    , m_velocity(control_interface->velocity_factor())
    , m_duration(control_interface->max_planning_duration())
    , m_iterations(control_interface->max_planning_iterations())
    , m_q(Eigen::Vector3f::Zero())
    , m_w({0.f, 0.f, 1.f})
    , m_eef_pos(Eigen::Vector3f::Zero())
    , m_eef_zyx(Eigen::Vector3f::Zero())
    , m_jog_pos_displacement(Eigen::Vector3f::Zero())
    , m_jog_zyx_displacement(Eigen::Vector3f::Zero())
    , m_trajectory(Trajectory::PTP, {std::make_pair(Trajectory::PTP, "P2P"), std::make_pair(Trajectory::LIN, "LIN"), std::make_pair(Trajectory::SCREW, "SCREW")})
    , m_control_mode(ControlMode::SIMULATION, {std::make_pair(ControlMode::PREVIEW, "PREVIEW"), std::make_pair(ControlMode::SIMULATION, "SIMULATION")})
    , m_control_space(ControlSpace::JOINT, {std::make_pair(ControlSpace::JOINT, "JOINT SPACE"), std::make_pair(ControlSpace::TASK, "TASK SPACE"), std::make_pair(ControlSpace::TOOL, "TOOL FRAME")})
    , m_control_widgets(ControlWidgets::FIELDS, {std::make_pair(ControlWidgets::SLIDERS, "SLIDERS"), std::make_pair(ControlWidgets::FIELDS, "FIELDS")})
    , m_control_interface(std::move(control_interface))
{
    m_joint_labels.resize(m_control_interface->joint_count());
    for(auto i = 0u; i < m_control_interface->joint_count(); i++)
        m_joint_labels[i] = "j" + std::to_string(i + 1);

    m_joint_pos = m_control_interface->joint_positionsf() * utility::rad_to_deg;
}

RobotControlWindows::RobotControlWindows(std::shared_ptr<Simulation::RobotControlInterface> control_interface, State state)
    : m_theta(0.f)
    , m_pitch(0.f)
    , m_velocity(state.velocity)
    , m_duration(state.duration)
    , m_iterations(state.iterations)
    , m_q(Eigen::Vector3f::Zero())
    , m_w({0.f, 0.f, 1.f})
    , m_eef_pos(Eigen::Vector3f::Zero())
    , m_eef_zyx(Eigen::Vector3f::Zero())
    , m_jog_pos_displacement(Eigen::Vector3f::Zero())
    , m_jog_zyx_displacement(Eigen::Vector3f::Zero())
    , m_trajectory(state.trajectory, {std::make_pair(Trajectory::PTP, "P2P"), std::make_pair(Trajectory::LIN, "LIN"), std::make_pair(Trajectory::SCREW, "SCREW")})
    , m_control_mode(state.control_mode, {std::make_pair(ControlMode::PREVIEW, "PREVIEW"), std::make_pair(ControlMode::SIMULATION, "SIMULATION")})
    , m_control_space(state.control_space, {std::make_pair(ControlSpace::JOINT, "JOINT SPACE"), std::make_pair(ControlSpace::TASK, "TASK SPACE"), std::make_pair(ControlSpace::TOOL, "TOOL FRAME")})
    , m_control_widgets(ControlWidgets::FIELDS, {std::make_pair(ControlWidgets::SLIDERS, "SLIDERS"), std::make_pair(ControlWidgets::FIELDS, "FIELDS")})
    , m_control_interface(std::move(control_interface))
{
    m_joint_labels.resize(m_control_interface->joint_count());
    for(auto i = 0u; i < m_control_interface->joint_count(); i++)
        m_joint_labels[i] = "j" + std::to_string(i + 1);

    m_joint_pos = m_control_interface->joint_positionsf() * utility::rad_to_deg;
}

RobotControlWindows::State RobotControlWindows::state() const
{
    return State{
        m_duration,
        m_iterations,
        m_velocity,
        m_trajectory.value(),
        m_control_mode.value(),
        m_control_space.value()
    };
}

void RobotControlWindows::render()
{
    ImGui::Begin("Robot controls");
    render_control_options();
    ImGui::End();

    ImGui::Begin("Control parameters");
    render_control_parameters();
    ImGui::End();

    ImGui::Begin("Control inputs");
    render_control_inputs();
    ImGui::End();
}

void RobotControlWindows::render_control_options()
{
    int selected_mode_item = m_control_mode.value_index();
    int selected_space_item = m_control_space.value_index();

    const auto &control_labels = m_control_mode.labels();
    if(ImGui::Combo("Control mode", &selected_mode_item, control_labels.data(), control_labels.size()))
    {
        m_control_mode.set(control_labels[selected_mode_item]);
    }

    const auto &space_labels = m_control_space.labels();
    if(ImGui::Combo("Control space", &selected_space_item, space_labels.data(), space_labels.size()))
    {
        m_control_space.set(space_labels[selected_space_item]);
        if(m_control_space != ControlSpace::JOINT && m_control_mode == ControlMode::PREVIEW)
            reset_to_current_pose();
    }

    if(m_control_space == ControlSpace::JOINT)
    {
        int widgets_selection = m_control_widgets.value_index();
        const auto &labels = m_control_widgets.labels();
        if(ImGui::Combo("Widgets", &widgets_selection, labels.data(), labels.size()))
            m_control_widgets.set(labels[widgets_selection]);
    }
    else if(m_control_space == ControlSpace::TASK)
    {
        int trajectory_selection = m_trajectory.value_index();
        const auto &labels = m_trajectory.labels();
        if(ImGui::Combo("Trajectory", &trajectory_selection, labels.data(), labels.size()))
            m_trajectory.set(labels[trajectory_selection]);
    }
}

void RobotControlWindows::render_control_parameters()
{
    if(ImGui::SliderInt("IK duration (msec)", &m_duration, 1, 1000))
        m_control_interface->max_planning_duration_changed(m_duration);
    if(ImGui::InputInt("IK duration (msec) ", &m_duration))
    {
        m_duration = std::abs(m_duration);
        m_control_interface->max_planning_duration_changed(m_duration);
    }
    ImGui::Separator();
    if(ImGui::SliderInt("IK iterations", &m_iterations, 1, 50000u))
        m_control_interface->max_planning_iterations_changed(m_iterations);
    if(ImGui::InputInt("IK iterations ", &m_iterations))
    {
        m_iterations = std::abs(m_iterations);
        m_control_interface->max_planning_iterations_changed(m_iterations);
    }
    ImGui::Separator();
    if(ImGui::SliderFloat("Velocity factor", &m_velocity, 0.f, 1.f))
        m_control_interface->velocity_factor_changed(m_velocity);
    auto limits = m_control_interface->joint_limits();
    for(auto i = 0u; i < m_control_interface->joint_count(); i++)
    {
        std::string idx = "J" + std::to_string(i + 1u);
        ImGui::Value((idx + " v_max").c_str(), static_cast<float>(limits.velocity[i] * utility::rad_to_deg) * m_velocity);
        ImGui::Value((idx + " a_max").c_str(), static_cast<float>(limits.acceleration[i] * utility::rad_to_deg) * m_velocity);
        if(i != m_control_interface->joint_count() - 1u)
            ImGui::Separator();
    }
}

void RobotControlWindows::render_control_inputs()
{
    if(m_control_space == ControlSpace::JOINT)
        render_joint_space();
    else if(m_control_mode == ControlMode::PREVIEW)
    {
        if(m_control_space == ControlSpace::TASK)
        {
            if(m_trajectory == Trajectory::SCREW)
                render_task_space_screw();
            else
                render_task_space_preview();
        }
        else if(m_control_space == ControlSpace::TOOL)
        {
            render_tool_frame_preview();
        }
    }
    else if(m_control_space == ControlSpace::TOOL)
    {
        render_tool_frame_preview();
    }
    else if(m_control_space == ControlSpace::TASK)
    {
        if(m_trajectory == Trajectory::SCREW)
            render_task_space_screw();
        else
            render_task_space_lin_p2p();
    }
}

void RobotControlWindows::render_joint_space()
{
    auto limits = m_control_interface->joint_limits();
    for(auto i = 0; i < m_control_interface->joint_count(); i++)
    {
        if(m_control_widgets == ControlWidgets::SLIDERS)
            ImGui::SliderFloat(m_joint_labels[i].c_str(), &m_joint_pos[i], limits.lower_position[i] * utility::rad_to_deg, limits.upper_position[i] * utility::rad_to_deg);
        else
            ImGui::InputFloat(m_joint_labels[i].c_str(), &m_joint_pos[i], limits.lower_position[i] * utility::rad_to_deg, limits.upper_position[i] * utility::rad_to_deg);
    }

    if(ImGui::Button("Reset to zero"))
        m_joint_pos = Eigen::VectorXf::Zero(6);
    if(m_control_mode == ControlMode::PREVIEW)
        m_control_interface->on_preview_joint_configuration(utility::to_eigen_vectord(m_joint_pos) * utility::deg_to_rad);
    else
    {
        if(ImGui::Button("Reset to current"))
            m_joint_pos = utility::to_eigen_vectorf(m_control_interface->joint_positions()) * utility::rad_to_deg;
        ImGui::SameLine();
        if(ImGui::Button("Move"))
            m_control_interface->on_joint_space_trajectory({utility::to_eigen_vectord(m_joint_pos) * utility::deg_to_rad});
    }
}

void RobotControlWindows::render_task_space_preview()
{
    auto preview_setpoint = [&]()
    {
        m_control_interface->on_preview_task_space_pose(to_matrix(m_eef_pos, m_eef_zyx * utility::deg_to_rad));
    };
    const char *position_labels[3] = {"X", "Y", "Z"};
    const char *orientation_labels[3] = {"A", "B", "C"};
    auto on_change = [&](int) { preview_setpoint(); };
    Widgets::render_float3_slider(m_eef_pos, position_labels, -1.f, 1.f, on_change);
    ImGui::NewLine();
    Widgets::render_float3_slider_with_reset(m_eef_zyx, orientation_labels, -180.f, 180.f, on_change);
}

void RobotControlWindows::render_tool_frame_preview()
{
    const char *position_labels[3] = {"X", "Y", "Z"};
    const char *orientation_labels[3] = {"A", "B", "C"};
    ImGui::InputFloat3("XYZ", &m_eef_pos[0]);
    ImGui::InputFloat3("ABC", &m_eef_zyx[0]);
    if(ImGui::Button("Reset"))
        reset_to_current_pose();
    auto on_change = [&](int)
    {
        if(m_control_mode == ControlMode::PREVIEW)
            m_control_interface->on_preview_tool_frame_jog(
                to_matrix(m_eef_pos, m_eef_zyx * utility::deg_to_rad),
                utility::to_eigen_vectord(m_jog_pos_displacement),
                utility::to_eigen_vectord(m_jog_zyx_displacement * utility::deg_to_rad)
            );
    };
    Widgets::render_float3_slider(m_jog_pos_displacement, position_labels, -1.f, 1.f, on_change);
    Widgets::render_float3_slider(m_jog_zyx_displacement, orientation_labels, -180.f, 180.f, on_change);
}

void RobotControlWindows::render_task_space_lin_p2p()
{
    const char *position_labels[3] = {"X", "Y", "Z"};
    const char *orientation_labels[3] = {"A", "B", "C"};
    auto no_change = [&](int)
    {
    };
    Widgets::render_float3_inputs(m_eef_pos, position_labels, 0.01f, 0.1f, no_change);
    ImGui::NewLine();
    Widgets::render_float3_inputs(m_eef_zyx, orientation_labels, 0.01f, 0.1f, no_change);
    if(ImGui::Button("Reset to current"))
        reset_to_current_pose();
    ImGui::SameLine();
    if(ImGui::Button("Move"))
    {
        if(m_trajectory == Trajectory::LIN)
            m_control_interface->on_task_space_lin(utility::to_eigen_vectord(m_eef_pos), utility::to_eigen_vectord(m_eef_zyx) * utility::deg_to_rad);
        else
            m_control_interface->on_task_space_ptp(utility::to_eigen_vectord(m_eef_pos), utility::to_eigen_vectord(m_eef_zyx) * utility::deg_to_rad);
    }
}

void RobotControlWindows::render_task_space_screw()
{
    const char *position_labels[3] = {"X", "Y", "Z"};
    const char *orientation_labels[3] = {"A", "B", "C"};
    ImGui::InputFloat3("XYZ", &m_eef_pos[0]);
    ImGui::InputFloat3("ABC", &m_eef_zyx[0]);
    if(ImGui::Button("Reset start"))
        reset_to_current_pose();

    ImGui::InputFloat3("q", m_q.data());
    ImGui::InputFloat3("w", m_w.data());
    ImGui::SameLine();
    if(ImGui::Button("Normalize"))
        m_w.normalize();
    ImGui::InputFloat("h", &m_pitch);
    if(ImGui::SliderFloat("theta", &m_theta, -180.f, 180.f))
    {
        if(m_control_mode == ControlMode::PREVIEW)
            m_control_interface->on_preview_task_space_screw(
                to_matrix(m_eef_pos, m_eef_zyx * utility::deg_to_rad), utility::to_eigen_vectord(m_w).normalized(), utility::to_eigen_vectord(m_q), m_theta * utility::deg_to_rad, m_pitch);
    }
    if(ImGui::Button("Reset"))
    {
        m_q = Eigen::Vector3f::Zero();
        m_w = {0.f, 0.f, 1.f};
        m_theta = 0.f;
        m_pitch = 0.f;
        m_control_interface->on_preview_task_space_screw(
            to_matrix(m_eef_pos, m_eef_zyx * utility::deg_to_rad), utility::to_eigen_vectord(m_w).normalized(), utility::to_eigen_vectord(m_q), m_theta * utility::deg_to_rad, m_pitch);
    }
    if(m_control_mode != ControlMode::PREVIEW)
    {
        ImGui::SameLine();
        if(ImGui::Button("Move"))
            m_control_interface->on_task_space_screw(utility::to_eigen_vectord(m_w).normalized(), utility::to_eigen_vectord(m_q), m_theta * utility::deg_to_rad, m_pitch);
    }
}

void RobotControlWindows::reset_to_current_pose()
{
    m_q = Eigen::Vector3f::Zero();
    m_w = Eigen::Vector3f{0.f, 0.f, 1.f};
    m_theta = 0.f;
    m_pitch = 0.f;
    m_jog_pos_displacement = Eigen::Vector3f::Zero();
    m_jog_zyx_displacement = Eigen::Vector3f::Zero();
    m_eef_pos = utility::to_eigen_vectorf(m_control_interface->eef_position());
    m_eef_zyx = utility::to_eigen_vectorf(m_control_interface->eef_orientation_zyx() * utility::rad_to_deg);
}

Eigen::Matrix4d RobotControlWindows::to_matrix(const Eigen::Vector3f &position, const Eigen::Vector3f &euler_zyx)
{
    return utility::transformation_matrix(
        utility::rotation_matrix_from_euler_zyx(utility::to_eigen_vectord(euler_zyx)),
        utility::to_eigen_vectord(position)
    );
}
