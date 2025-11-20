#include "visualization/imgui/imguiwindowcontext.h"

using namespace AIS4104::Visualization;

ImguiWindowContext::ImguiWindowContext(const threepp::Canvas &canvas, std::vector<std::shared_ptr<ImguiWindow>> windows)
    : ImguiFunctionalContext(canvas.windowPtr(), [this] { render_windows(); })
    , m_windows(std::move(windows))
{
}

void ImguiWindowContext::render_windows()
{
    for(auto &item : m_windows)
        item->render();
}

void ImguiWindowContext::initialize()
{
    for(auto &item : m_windows)
        item->initialize();
}
