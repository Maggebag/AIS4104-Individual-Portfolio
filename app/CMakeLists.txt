add_executable(app
        src/app/main.cpp
        src/app/robotwrapper.h
        src/app/robotwrapper.cpp
        src/app/hardcoded_ui_states.h

        src/app/urdf_loaders.h
        src/app/hardcoded_loaders.h

        src/app/implementations/motionplannerimpl.h
        src/app/implementations/motionplannerimpl.cpp
        src/app/implementations/screwskinematicssolver.h
        src/app/implementations/screwskinematicssolver.cpp
        src/app/implementations/ptptrajectorygenerator.hpp
        src/app/implementations/ptptrajectorygenerator.cpp

        src/app/examples/tracikkinematicssolver.h
        src/app/examples/tracikkinematicssolver.cpp
        src/app/examples/ptplinetrajectorygeneratorexample.h
        src/app/examples/ptplinetrajectorygeneratorexample.cpp
)

target_pedantic_build(app)
target_link_libraries(app PRIVATE visualization utility threepp)
target_include_directories(app PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)