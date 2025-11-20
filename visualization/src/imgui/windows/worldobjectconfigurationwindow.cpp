#include "visualization/imgui/windows/worldobjectconfigurationwindow.h"

#include "visualization/imgui/widgets.h"

#include <utility/math.h>

#include <utility/vectors.h>

#include <string.h>

using namespace AIS4104::Visualization;

WorldObjectConfigurationWindow::WorldObjectConfigurationWindow(RobotScene &scene)
    : m_active(false)
    , m_model_path('\0')
    , m_scene(scene)
    , m_gfx_scale{1.f, 1.f, 1.f}
    , m_gfx_offset{0.f, 0.f, 0.f}
    , m_gfx_euler_zyx{0.f, 0.f, 0.f}
    , m_world_view(WorldView::LOAD_STL, {std::make_pair(WorldView::TRANSFORM, "Transform"), std::make_pair(WorldView::LOAD_STL, "Load .stl")})
{
}

WorldObjectConfigurationWindow::WorldObjectConfigurationWindow(RobotScene &scene, State state)
    : m_active(state.active)
    , m_model_path('\0')
    , m_scene(scene)
    , m_gfx_scale{state.gfx_scale}
    , m_gfx_offset{state.gfx_offset}
    , m_gfx_euler_zyx{state.gfx_euler_zyx}
    , m_world_view(state.selected_view, {std::make_pair(WorldView::TRANSFORM, "Transform"), std::make_pair(WorldView::LOAD_STL, "Load .stl")})
{
    std::strncpy(m_model_path, &state.model_path[0], sizeof(m_model_path));
}

WorldObjectConfigurationWindow::State WorldObjectConfigurationWindow::state() const
{
    return State{
        m_active,
        m_model_path,
        m_world_view.value(),
        m_gfx_scale,
        m_gfx_offset,
        m_gfx_euler_zyx,
    };
}

void WorldObjectConfigurationWindow::render()
{
    ImGui::Begin("World object settings");
    if(m_world_object != nullptr)
    {
        if(ImGui::Checkbox("Active", &m_active))
        {
            if(m_active)
            {
                activate_loaded_world_object();
                m_world_view = WorldView::TRANSFORM;
            }
            else
            {
                deactivate_loaded_world_object();
                m_world_view = WorldView::LOAD_STL;
            }
        }
        ImGui::Text("Loaded model: %s", m_model_path);
    }

    if(m_world_view == WorldView::TRANSFORM)
        render_graphics_transform();
    else
        render_stl_loader();
    ImGui::End();
}

void WorldObjectConfigurationWindow::initialize()
{
    if(load_stl())
    {
        if(m_active)
        {
            activate_loaded_world_object();
            m_world_view = WorldView::TRANSFORM;
        }
    }
    else
        m_world_view = WorldView::LOAD_STL;
}

void WorldObjectConfigurationWindow::render_stl_loader()
{
    ImGui::InputText("STL file", m_model_path, 1024);
    if(ImGui::Button("Load"))
    {
        m_active = load_stl();
        if(m_active)
        {
            m_world_view = WorldView::TRANSFORM;
            activate_loaded_world_object();
        }
    }
}

void WorldObjectConfigurationWindow::render_graphics_transform()
{
    const char *scale_labels[3] = {"SX", "SY", "SZ"};
    const char *position_labels[3] = {"X", "Y", "Z"};
    const char *orientation_labels[3] = {"A", "B", "C"};
    Widgets::render_float3_inputs(m_gfx_offset, position_labels, 0.1f, 1.f, [&](int) { assign_gfx_transform(); });
    ImGui::NewLine();
    Widgets::render_float3_inputs_with_reset(m_gfx_euler_zyx, orientation_labels, 1.f, 10.f, [&](int) { assign_gfx_transform(); });
    ImGui::NewLine();
    Widgets::render_float3_inputs(m_gfx_scale, scale_labels, 0.1f, 1.f, [&](int) { assign_gfx_transform(); });
}

void WorldObjectConfigurationWindow::assign_gfx_transform()
{
    Eigen::Matrix4d tf = utility::transformation_matrix(utility::rotation_matrix_from_euler_zyx(utility::to_eigen_vectord(m_gfx_euler_zyx * utility::deg_to_rad)), utility::to_eigen_vectord(m_gfx_offset));
    threepp::Matrix4 world_object_tf;
    world_object_tf.set(
        static_cast<float>(tf(0, 0)), static_cast<float>(tf(0, 1)), static_cast<float>(tf(0, 2)), static_cast<float>(tf(0, 3)),
        static_cast<float>(tf(1, 0)), static_cast<float>(tf(1, 1)), static_cast<float>(tf(1, 2)), static_cast<float>(tf(1, 3)),
        static_cast<float>(tf(2, 0)), static_cast<float>(tf(2, 1)), static_cast<float>(tf(2, 2)), static_cast<float>(tf(2, 3)),
        static_cast<float>(tf(3, 0)), static_cast<float>(tf(3, 1)), static_cast<float>(tf(3, 2)), static_cast<float>(tf(3, 3))
    );
    threepp::Vector3 scale(m_gfx_scale.x(), m_gfx_scale.y(), m_gfx_scale.z());
    threepp::Vector3 offset(m_gfx_offset.x(), m_gfx_offset.y(), m_gfx_offset.z());
    m_world_object->scale = scale;
    m_world_object->position = offset;
    m_world_object->setRotationFromMatrix(world_object_tf);
}

void WorldObjectConfigurationWindow::activate_loaded_world_object()
{
    m_scene.set_world_object(m_world_object);
    assign_gfx_transform();
}

void WorldObjectConfigurationWindow::deactivate_loaded_world_object()
{
    m_scene.clear_world_object();
}

void WorldObjectConfigurationWindow::clear_world_object()
{
    deactivate_loaded_world_object();
    if(m_world_object)
        m_world_object.reset();
}

bool WorldObjectConfigurationWindow::load_stl()
{
    clear_world_object();
    threepp::STLLoader loader;
    std::string path(m_model_path);
    if(path.empty())
        return false;
    auto geometry = loader.load(path);
    if(geometry == nullptr)
        return false;
    auto material = threepp::MeshPhongMaterial::create({{"flatShading", true}, {"color", threepp::Color::gray}});
    auto mesh = threepp::Mesh::create(geometry, material);
    m_world_object = mesh;
    return true;
}
