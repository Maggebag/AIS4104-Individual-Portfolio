#include <utility>

#include "visualization/robotscene.h"

using namespace AIS4104::Visualization;

constexpr float frustum_size = 5.f;

RobotScene::RobotScene(std::shared_ptr<threepp::Robot> robot, std::vector<std::shared_ptr<ImguiWindow>> imgui_items, Projection projection)
    : m_projection(projection)
    , m_robot(std::move(robot))
{
    std::unordered_map<std::string, threepp::Canvas::ParameterValue> params{{"aa", 4}};
    m_canvas = std::make_unique<threepp::Canvas>("Robot renderer", params);
    m_imgui = std::make_shared<ImguiWindowContext>(*m_canvas, std::move(imgui_items));
    m_imgui->setFontScale(1.f);

    m_renderer = std::make_unique<threepp::GLRenderer>(m_canvas->size());
    m_renderer->setClearColor(threepp::Color::aliceblue);

    setup_scene();

    m_capture = std::make_unique<threepp::IOCapture>();
    m_capture->preventMouseEvent = []
    {
        return ImGui::GetIO().WantCaptureMouse;
    };
    m_canvas->setIOCapture(m_capture.get());
    m_canvas->onWindowResize([&](threepp::WindowSize size)
        {
            if(m_projection == Projection::PERSPECTIVE)
            {
                static_cast<threepp::PerspectiveCamera&>(*m_camera).aspect = size.aspect();
            }
            else
            {
                auto &camera = static_cast<threepp::OrthographicCamera&>(*m_camera);
                camera.left = -frustum_size * size.aspect() / 2.f;
                camera.right = frustum_size * size.aspect() / 2.f;
                camera.top = frustum_size / 2.f;
                camera.bottom = -frustum_size / 2.f;
            }
            m_camera->updateProjectionMatrix();
            m_renderer->setSize(size);
        }
    );
}

RobotScene::~RobotScene()
{
    if(m_canvas->isOpen())
        m_canvas->close();
    m_imgui.reset();
}

void RobotScene::show_grid(bool show)
{
    m_gridhelper->visible = show;
}

void RobotScene::set_background_color(threepp::Color color)
{
    m_scene->background = color;
}

void RobotScene::set_tool(std::shared_ptr<threepp::Object3D> tool)
{
    threepp::Matrix4 i;
    set_tool(std::move(tool), i.identity());
}

void RobotScene::set_tool(std::shared_ptr<threepp::Object3D> tool, threepp::Matrix4 tool_tf)
{
    reset_tool();
    m_tool = std::move(tool);
    m_scene->add(m_tool);
    m_tool_tf = tool_tf;
}

void RobotScene::set_tool_transform(threepp::Matrix4 tool_tf)
{
    m_tool_tf = tool_tf;
}

void RobotScene::reset_tool()
{
    if(m_tool)
    {
        m_scene->remove(*m_tool);
        m_tool.reset();
        m_tool_tf.identity();
    }
}

void RobotScene::set_world_object(std::shared_ptr<threepp::Object3D> world_object)
{
    clear_world_object();
    m_world_object = std::move(world_object);
    m_scene->add(m_world_object);
}

void RobotScene::clear_world_object()
{
    if(m_world_object)
    {
        m_scene->remove(*m_world_object);
        m_world_object.reset();
    }
}

threepp::Robot& RobotScene::robot()
{
    return *m_robot;
}

const threepp::Robot& RobotScene::robot() const
{
    return *m_robot;
}

void RobotScene::render() const
{
    m_imgui->initialize();
    threepp::Clock clock;
    m_canvas->animate([&]()
        {
            if(m_tool != nullptr)
            {
                threepp::Matrix4 t_sb = m_robot->computeEndEffectorTransform(m_robot->jointValuesWithConversionFromRadiansToDeg(), true);
                m_tool->position.setFromMatrixPosition(t_sb);
                m_tool->quaternion.setFromRotationMatrix(t_sb.multiply(m_tool_tf));
            }
            m_renderer->render(*m_scene, *m_camera);
            m_imgui->render();
        }
    );
}

void RobotScene::setup_scene()
{
    m_scene = threepp::Scene::create();
    add_scene_grid();
    add_camera();
    add_scene_objects();
}

void RobotScene::add_scene_grid()
{
    auto size = 10.f;
    auto material = threepp::ShadowMaterial::create();
    auto plane = threepp::Mesh::create(threepp::PlaneGeometry::create(size, size), material);
    plane->rotation.y = -threepp::math::PI / 2.f;
    plane->receiveShadow = true;

    m_gridhelper = threepp::GridHelper::create(size, size, threepp::Color::yellowgreen);
    m_gridhelper->rotation.y = threepp::math::PI / 2.f;
    plane->add(m_gridhelper);
    m_scene->add(m_gridhelper);
}

void RobotScene::add_camera()
{
    if(m_projection == Projection::PERSPECTIVE)
        m_camera = threepp::PerspectiveCamera::create(75.f, m_canvas->aspect(), 0.1f, 100.f);
    else
    {
        m_camera = threepp::OrthographicCamera::create(
            -(frustum_size * m_canvas->aspect()) / 2.f,
            (frustum_size * m_canvas->aspect()) / 2.f,
            frustum_size / 2.f,
            -frustum_size / 2.f,
            0.001f,
            100.f
        );
    }
    m_camera->position.set(0.f, 1.25f, 1.f);
    m_controls = std::make_unique<threepp::OrbitControls>(*m_camera, *m_canvas);
}

void RobotScene::add_scene_objects()
{
    auto light = threepp::HemisphereLight::create(threepp::Color::aliceblue, threepp::Color::grey);
    m_scene->add(light);
    m_scene->add(*m_robot);
    // m_scene->add(*m_tool);
}
