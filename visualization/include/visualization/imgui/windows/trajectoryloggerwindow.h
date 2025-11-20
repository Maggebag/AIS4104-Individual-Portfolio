#ifndef AIS4104_VISUALIZATION_TRAJECTORYLOGGERWINDOW_H
#define AIS4104_VISUALIZATION_TRAJECTORYLOGGERWINDOW_H

#include <simulation/robotcontrolinterface.h>

#include <visualization/imgui/imguiwindowcontext.h>

namespace AIS4104::Visualization {

class TrajectoryLoggerWindow : public ImguiWindow
{
public:
    typedef Simulation::LoggerParameters State;

    TrajectoryLoggerWindow(Simulation::RobotControlInterface &control_interface);
    TrajectoryLoggerWindow(Simulation::RobotControlInterface &control_interface, State state);

    State state() const;

    void render() override;

private:
    bool m_active;
    char m_directory_path[1024];
    Simulation::RobotControlInterface &m_control_interface;

    void update_folder();
};

}

#endif
