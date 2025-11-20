#include "visualization/imgui/windows/trajectoryloggerwindow.h"

#include "spdlog/spdlog.h"

using namespace AIS4104::Visualization;

TrajectoryLoggerWindow::TrajectoryLoggerWindow(Simulation::RobotControlInterface &control_interface)
    : TrajectoryLoggerWindow(control_interface, control_interface.logger_parameters())
{
}

TrajectoryLoggerWindow::TrajectoryLoggerWindow(Simulation::RobotControlInterface &control_interface, State state)
    : m_active(false)
    , m_directory_path('\0')
    , m_control_interface(control_interface)
{
    m_active = state.active;
    std::string directory = state.directory.string();
    std::strncpy(m_directory_path, &directory[0], sizeof(m_directory_path));
    m_control_interface.set_logger_parameters(state);
}

TrajectoryLoggerWindow::State TrajectoryLoggerWindow::state() const
{
    return State{m_active, m_directory_path};
}

void TrajectoryLoggerWindow::render()
{
    ImGui::Begin("Logging");
    if(ImGui::Checkbox("Active", &m_active))
        m_control_interface.set_logger_parameters(state());
    ImGui::InputText("Output folder", m_directory_path, 1024);
    if(ImGui::Button("Update folder"))
        update_folder();
    ImGui::End();
}

void TrajectoryLoggerWindow::update_folder()
{
    std::filesystem::path directory(m_directory_path);
    if(!exists(directory))
        if(!create_directories(directory))
        {
            spdlog::error("Unable to create directory: {}", directory.string());
            return;
        }
    m_control_interface.set_logger_parameters(state());
}
