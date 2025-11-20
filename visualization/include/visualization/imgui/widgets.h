#ifndef AIS4104_VISUALIZATION_IMGUI_WIDGETS_H
#define AIS4104_VISUALIZATION_IMGUI_WIDGETS_H

#include "visualization/imgui/imguiwindowcontext.h"

#include <Eigen/Dense>

namespace AIS4104::Visualization::Widgets {

static void render_float3_inputs(Eigen::Vector3f &vector, const char *labels[3], float step, float step_fast, const std::function<void(int idx)> &on_change)
{
    if(ImGui::InputFloat(labels[0], &vector[0], step, step_fast))
        on_change(0);
    if(ImGui::InputFloat(labels[1], &vector[1], step, step_fast))
        on_change(1);
    if(ImGui::InputFloat(labels[2], &vector[2], step, step_fast))
        on_change(2);
}

static void render_float_inputs_with_reset(float &f, const char *label, float step, float step_fast, const std::function<void()> &on_change)
{
    if(ImGui::InputFloat(label, &f, step, step_fast))
        on_change();
    ImGui::SameLine();
    if(ImGui::Button((std::string("Reset ") + label).c_str()))
    {
        f = 0.f;
        on_change();
    }
}

static void render_float3_inputs_with_reset(Eigen::Vector3f &vector, const char *labels[3], float step, float step_fast, const std::function<void(int idx)> &on_change)
{
    render_float_inputs_with_reset(vector[0], labels[0], step, step_fast, [&]() { on_change(0); });
    render_float_inputs_with_reset(vector[1], labels[1], step, step_fast, [&]() { on_change(1); });
    render_float_inputs_with_reset(vector[2], labels[2], step, step_fast, [&]() { on_change(2); });
}

static void render_float_slider_with_reset(float &f, const char *label, float min, float max, const std::function<void()> &on_change)
{
    if(ImGui::SliderFloat(label, &f, min, max))
        on_change();
    ImGui::SameLine();
    if(ImGui::Button((std::string("Reset ") + label).c_str()))
    {
        f = 0.f;
        on_change();
    }
}

static void render_float3_slider(Eigen::Vector3f &vector, const char *labels[3], float min, float max, const std::function<void(int idx)> &on_change)
{
    if(ImGui::SliderFloat(labels[0], &vector[0], min, max))
        on_change(0);
    if(ImGui::SliderFloat(labels[1], &vector[1], min, max))
        on_change(1);
    if(ImGui::SliderFloat(labels[2], &vector[2], min, max))
        on_change(2);
}

static void render_float3_slider_with_reset(Eigen::Vector3f &vector, const char *labels[3], float min, float max, const std::function<void(int idx)> &on_change)
{
    render_float_slider_with_reset(vector[0], labels[0], min, max, [&]() { on_change(0); });
    render_float_slider_with_reset(vector[1], labels[1], min, max, [&]() { on_change(1); });
    render_float_slider_with_reset(vector[2], labels[2], min, max, [&]() { on_change(2); });
}

}

#endif
