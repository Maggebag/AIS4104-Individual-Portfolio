#include "visualization/imgui/windows/toolwindow.h"

#include <utility/math.h>
#include <utility/vectors.h>

#include <visualization/imgui/widgets.h>

using namespace AIS4104::Visualization;

constexpr double rad2deg = 57.29577951308232;
constexpr double deg2rad = 0.0174532925199433;

ToolConfigurationWindow::ToolConfigurationWindow(RobotScene &scene, std::shared_ptr<Simulation::Robot> robot)
    : m_active(false)
    , m_model_path('\0')
    , m_gfx_scale{1.f, 1.f, 1.f}
    , m_gfx_offset{0.f, 0.f, 0.f}
    , m_gfx_euler_zyx{0.f, 0.f, 0.f}
    , m_kinematic_offset{0.0f, 0.f, 0.f}
    , m_kinematic_euler_zyx{0.f, 0.f, 0.f}
    , m_scene(scene)
    , m_tool_view(ToolView::LOAD_STL, {std::make_pair(ToolView::KINEMATICS_TFX, "Kinematics transform"), std::make_pair(ToolView::GFX_TFX, "Graphics transform"), std::make_pair(ToolView::LOAD_STL, "Load .stl")})
    , m_robot(std::move(robot))
{
}

ToolConfigurationWindow::ToolConfigurationWindow(RobotScene &scene, std::shared_ptr<Simulation::Robot> robot, State state)
    : m_active(state.active)
    , m_model_path('\0')
    , m_gfx_scale(state.gfx_scale)
    , m_gfx_offset(state.gfx_offset)
    , m_gfx_euler_zyx(state.gfx_euler_zyx)
    , m_kinematic_offset(state.kinematic_offset)
    , m_kinematic_euler_zyx(state.kinematic_euler_zyx)
    , m_scene(scene)
    , m_tool_view(state.selected_view, {std::make_pair(ToolView::KINEMATICS_TFX, "Kinematics transform"), std::make_pair(ToolView::GFX_TFX, "Graphics transform"), std::make_pair(ToolView::LOAD_STL, "Load .stl")})
    , m_robot(std::move(robot))
{
    std::strncpy(m_model_path, &state.model_path[0], sizeof(m_model_path));
}

ToolConfigurationWindow::State ToolConfigurationWindow::state() const
{
    return State{
        m_active,
        m_model_path,
        m_tool_view.value(),
        m_gfx_scale,
        m_gfx_offset,
        m_gfx_euler_zyx,
        m_kinematic_offset,
        m_kinematic_euler_zyx
    };
}

void ToolConfigurationWindow::render()
{
    ImGui::Begin("Tool settings");
    if(m_tool != nullptr)
    {
        if(ImGui::Checkbox("Active", &m_active))
        {
            if(!m_active)
                activate_default_tool();
            else
                activate_custom_tool();
        }
        ImGui::Text("Loaded model: %s", m_model_path);
    }

    int selected_type = m_tool_view.value_index();
    auto type_labels = m_tool_view.labels();
    if(ImGui::Combo("Tool view", &selected_type, type_labels.data(), type_labels.size()))
        m_tool_view.set(type_labels[selected_type]);
    if(m_tool_view == ToolView::GFX_TFX)
        render_graphics_transform();
    else if(m_tool_view == ToolView::KINEMATICS_TFX)
        render_kinematics_transform();
    else if(m_tool_view == ToolView::LOAD_STL)
        render_stl_loader();
    ImGui::End();
}

void ToolConfigurationWindow::initialize()
{
    if(load_stl() && m_active)
    {
        activate_custom_tool();
        m_tool_view = ToolView::KINEMATICS_TFX;
    }
    else
    {
        activate_default_tool();
        m_tool_view = ToolView::LOAD_STL;
    }
}

void ToolConfigurationWindow::render_stl_loader()
{
    ImGui::InputText("STL file", m_model_path, 1024);
    if(ImGui::Button("Load"))
    {
        m_active = load_stl();
        if(m_active)
        {
            m_tool_view = ToolView::KINEMATICS_TFX;
            activate_custom_tool();
        }
    }
}

void ToolConfigurationWindow::render_kinematics_transform()
{
    const char *position_labels[3] = {"X", "Y", "Z"};
    const char *orientation_labels[3] = {"A", "B", "C"};
    auto on_change = [&](int)
    {
    };
    Widgets::render_float3_inputs(m_kinematic_offset, position_labels, 0.01f, 0.1f, on_change);
    ImGui::NewLine();
    Widgets::render_float3_inputs_with_reset(m_kinematic_euler_zyx, orientation_labels, 1.f, 10.f, on_change);
    if(ImGui::Button("Set"))
        assign_kinematics_transform();
}

void ToolConfigurationWindow::render_graphics_transform()
{
    const char *scale_labels[3] = {"SX", "SY", "SZ"};
    const char *position_labels[3] = {"X", "Y", "Z"};
    const char *orientation_labels[3] = {"A", "B", "C"};
    auto on_change = [&](int)
    {
    };
    Widgets::render_float3_inputs(m_gfx_offset, position_labels, 0.01f, 0.1f, on_change);
    ImGui::NewLine();
    Widgets::render_float3_inputs_with_reset(m_gfx_euler_zyx, orientation_labels, 0.01f, 0.1f, on_change);
    ImGui::NewLine();
    Widgets::render_float3_inputs(m_gfx_scale, scale_labels, 0.01f, 0.1f, on_change);
    if(ImGui::Button("Set"))
        assign_gfx_transform();
}

void ToolConfigurationWindow::assign_gfx_transform()
{
    Eigen::Matrix4d tf = utility::transformation_matrix(utility::rotation_matrix_from_euler_zyx(utility::to_eigen_vectord(m_gfx_euler_zyx * deg2rad)), utility::to_eigen_vectord(m_gfx_offset));
    threepp::Matrix4 tool_tf;
    tool_tf.set(
        static_cast<float>(tf(0, 0)), static_cast<float>(tf(0, 1)), static_cast<float>(tf(0, 2)), static_cast<float>(tf(0, 3)),
        static_cast<float>(tf(1, 0)), static_cast<float>(tf(1, 1)), static_cast<float>(tf(1, 2)), static_cast<float>(tf(1, 3)),
        static_cast<float>(tf(2, 0)), static_cast<float>(tf(2, 1)), static_cast<float>(tf(2, 2)), static_cast<float>(tf(2, 3)),
        static_cast<float>(tf(3, 0)), static_cast<float>(tf(3, 1)), static_cast<float>(tf(3, 2)), static_cast<float>(tf(3, 3))
    );
    threepp::Vector3 scale(m_gfx_scale.x(), m_gfx_scale.y(), m_gfx_scale.z());
    m_tool->scale = scale;
    m_scene.set_tool_transform(tool_tf);
}

void ToolConfigurationWindow::assign_kinematics_transform()
{
    m_robot->set_tool_transform(utility::transformation_matrix(utility::rotation_matrix_from_euler_zyx(utility::to_eigen_vectord(m_kinematic_euler_zyx * deg2rad)), utility::to_eigen_vectord(m_kinematic_offset)));
}

void ToolConfigurationWindow::clear_assigned_kinematics_transform()
{
    m_robot->set_tool_transform(Eigen::Matrix4d::Identity());
}

void ToolConfigurationWindow::clear_tool()
{
    m_scene.reset_tool();
    if(m_tool)
        m_tool.reset();
}

void ToolConfigurationWindow::activate_custom_tool()
{
    m_scene.set_tool(m_tool);
    assign_gfx_transform();
    assign_kinematics_transform();
}

void ToolConfigurationWindow::activate_default_tool()
{
    threepp::Box3 robot_bounding_box;
    robot_bounding_box.setFromObject(m_scene.robot());
    auto axes = threepp::AxesHelper::create(robot_bounding_box.getSize().length() * 0.1f);
    m_scene.set_tool(axes, threepp::Matrix4().makeRotationY(-threepp::math::PI / 2.f));
    m_robot->set_tool_transform(Eigen::Matrix4d::Identity());
}

bool ToolConfigurationWindow::load_stl()
{
    threepp::STLLoader loader;
    std::string path(m_model_path);
    if(path.empty())
        return false;
    auto geometry = loader.load(path);
    if(geometry == nullptr)
        return false;
    auto material = threepp::MeshPhongMaterial::create({{"flatShading", true}, {"color", threepp::Color::gray}});
    auto mesh = threepp::Mesh::create(geometry, material);
    m_tool = mesh;
    return true;
}
