#ifndef AIS4104_VISUALIZATION_POSEDISPLAYWINDOW_H
#define AIS4104_VISUALIZATION_POSEDISPLAYWINDOW_H

#include <simulation/robotcontrolinterface.h>

#include <visualization/imgui/imguiwindowcontext.h>

#include <utility/enum.h>

namespace AIS4104::Visualization {

class PoseDisplayWindow : public ImguiWindow
{
    enum class FrameView { FLANGE, TOOL };

public:
    PoseDisplayWindow(const Simulation::RobotControlInterface &control_interface);

    void render() override;

private:
    utility::Enum<FrameView, 2> m_frame_view;
    const Simulation::RobotControlInterface &m_control_interface;
};

}

#endif
