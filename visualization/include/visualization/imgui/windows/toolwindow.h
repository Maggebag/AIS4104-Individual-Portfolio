#ifndef AIS4104_VISUALIZATION_TOOLWINDOW_H
#define AIS4104_VISUALIZATION_TOOLWINDOW_H

#include <simulation/robot.h>

#include <utility/enum.h>

#include <visualization/robotscene.h>

#include <visualization/imgui/imguiwindowcontext.h>

namespace AIS4104::Visualization {

class ToolConfigurationWindow : public ImguiWindow
{
public:
    enum class ToolView { KINEMATICS_TFX, GFX_TFX, LOAD_STL };

    struct State
    {
        bool active;
        std::string model_path;
        ToolView selected_view;
        Eigen::Vector3f gfx_scale;
        Eigen::Vector3f gfx_offset;
        Eigen::Vector3f gfx_euler_zyx;
        Eigen::Vector3f kinematic_offset;
        Eigen::Vector3f kinematic_euler_zyx;
    };

    ToolConfigurationWindow(RobotScene &scene, std::shared_ptr<Simulation::Robot> robot);
    ToolConfigurationWindow(RobotScene &scene, std::shared_ptr<Simulation::Robot> robot, State state);

    State state() const;

    void render() override;

    void initialize() override;

private:
    bool m_active;
    char m_model_path[1024];
    Eigen::Vector3f m_gfx_scale;
    Eigen::Vector3f m_gfx_offset;
    Eigen::Vector3f m_gfx_euler_zyx;
    Eigen::Vector3f m_kinematic_offset;
    Eigen::Vector3f m_kinematic_euler_zyx;
    RobotScene &m_scene;
    utility::Enum<ToolView, 3> m_tool_view;
    std::shared_ptr<threepp::Object3D> m_tool;
    std::shared_ptr<Simulation::Robot> m_robot;

    void render_stl_loader();
    void render_graphics_transform();
    void render_kinematics_transform();

    void assign_gfx_transform();
    void assign_kinematics_transform();
    void clear_assigned_kinematics_transform();

    void clear_tool();
    void activate_custom_tool();
    void activate_default_tool();

    bool load_stl();
};

}

#endif
