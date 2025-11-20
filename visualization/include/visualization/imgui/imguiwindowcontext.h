#ifndef AIS4104_IMGUIWINDOWCONTEXT_H
#define AIS4104_IMGUIWINDOWCONTEXT_H

#include <threepp/canvas/Canvas.hpp>

#include <threepp/extras/imgui/ImguiContext.hpp>

#include <memory>

namespace AIS4104::Visualization {

class ImguiWindow
{
public:
    ImguiWindow() = default;
    virtual ~ImguiWindow() = default;

    virtual void render() = 0;

    virtual void initialize()
    {
    }
};

class ImguiWindowContext : public ImguiFunctionalContext
{
public:
    ImguiWindowContext(const threepp::Canvas &canvas, std::vector<std::shared_ptr<ImguiWindow>> windows);

    void render_windows();
    void initialize();

private:
    std::vector<std::shared_ptr<ImguiWindow>> m_windows;
};

}

#endif
