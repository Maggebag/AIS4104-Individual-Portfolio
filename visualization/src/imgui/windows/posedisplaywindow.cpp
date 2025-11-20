#include "visualization/imgui/windows/posedisplaywindow.h"

#include <utility/math.h>
#include <utility/vectors.h>

using namespace AIS4104::Visualization;

namespace tpp = threepp;

PoseDisplayWindow::PoseDisplayWindow(const Simulation::RobotControlInterface &control_interface)
    : m_frame_view(FrameView::TOOL, {std::make_pair(FrameView::TOOL, "Tool"), std::make_pair(FrameView::FLANGE, "Flange")})
    , m_control_interface(control_interface)
{
}

void PoseDisplayWindow::render()
{
    ImGui::Begin("Pose");
    int selected_type = m_frame_view.value_index();
    auto type_labels = m_frame_view.labels();
    if(ImGui::Combo("Frame", &selected_type, type_labels.data(), type_labels.size()))
    {
        m_frame_view.set(type_labels[selected_type]);
    }
    Eigen::Matrix4d tf = Eigen::Matrix4d::Identity();
    if(m_frame_view == FrameView::FLANGE)
        tf = m_control_interface.flange_pose();
    else
        tf = m_control_interface.eef_pose();

    Eigen::Vector3f pos = utility::to_eigen_vectorf(tf.block<3, 1>(0, 3));
    Eigen::Vector3f zyx = utility::to_eigen_vectorf(utility::euler_zyx_from_rotation_matrix(tf.block<3, 3>(0, 0)) * utility::rad_to_deg);
    ImGui::Value("X", pos.x());
    ImGui::SameLine();
    ImGui::Value("A", zyx.x());

    ImGui::Value("Y", pos.y());
    ImGui::SameLine();
    ImGui::Value("B", zyx.y());

    ImGui::Value("Z", pos.z());
    ImGui::SameLine();
    ImGui::Value("C", zyx.z());
    ImGui::End();
}
