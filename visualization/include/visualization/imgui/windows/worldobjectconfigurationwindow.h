#ifndef AIS4104_VISUALIZATION_WORLDOBJECTCONFIGURATIONWINDOW_H
#define AIS4104_VISUALIZATION_WORLDOBJECTCONFIGURATIONWINDOW_H

#include "visualization/robotscene.h"

#include <utility/enum.h>

#include <Eigen/Dense>

namespace AIS4104::Visualization {

class WorldObjectConfigurationWindow : public ImguiWindow
{
public:
    enum class WorldView { LOAD_STL, TRANSFORM };

    struct State
    {
        bool active;
        std::string model_path;
        WorldView selected_view;
        Eigen::Vector3f gfx_scale = Eigen::Vector3f::Ones();
        Eigen::Vector3f gfx_offset = Eigen::Vector3f::Zero();
        Eigen::Vector3f gfx_euler_zyx = Eigen::Vector3f::Zero();
    };

    WorldObjectConfigurationWindow(RobotScene &scene);
    WorldObjectConfigurationWindow(RobotScene &scene, State state);

    State state() const;

    void render() override;

    void initialize() override;

private:
    bool m_active;
    char m_model_path[1024];
    RobotScene &m_scene;
    Eigen::Vector3f m_gfx_scale;
    Eigen::Vector3f m_gfx_offset;
    Eigen::Vector3f m_gfx_euler_zyx;
    utility::Enum<WorldView, 2> m_world_view;
    std::shared_ptr<threepp::Object3D> m_world_object;

    void render_stl_loader();

    void render_graphics_transform();

    void assign_gfx_transform();

    void activate_loaded_world_object();
    void deactivate_loaded_world_object();
    void clear_world_object();

    bool load_stl();
};

}

#endif
