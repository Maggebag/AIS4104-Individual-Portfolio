find_package(glfw3 CONFIG REQUIRED)

add_library(imgui STATIC)
add_library(imgui::imgui ALIAS imgui)

target_sources(imgui PRIVATE
        "src/imgui.cpp"
        "src/imgui_draw.cpp"
        "src/imgui_tables.cpp"
        "src/imgui_widgets.cpp"

        "src/backends/imgui_impl_glfw.cpp"
        "src/backends/imgui_impl_opengl3.cpp"
)

target_link_libraries(imgui PRIVATE glfw)
target_include_directories(imgui PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}/src/backends")
