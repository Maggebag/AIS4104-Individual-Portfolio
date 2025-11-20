find_package(Eigen3 CONFIG REQUIRED)

set(PUBLIC_HEADERS
    include/simulation/robot.h
    include/simulation/factory.h
    include/simulation/motionplanner.h
    include/simulation/kinematicssolver.h
    include/simulation/trajectorygenerator.h
    include/simulation/robotcontrolinterface.h
)

add_library(simulation
    src/simulation/factory.cpp
    src/simulation/kinematicssolver.cpp
    src/simulation/trajectoryexecutor.h
    src/simulation/trajectoryexecutor.cpp
    src/simulation/robotcontroldelegator.h
    src/simulation/robotcontroldelegator.cpp

    ${PUBLIC_HEADERS}
)

target_pedantic_build(simulation)
target_link_libraries(simulation PUBLIC utility)
target_include_directories(simulation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
