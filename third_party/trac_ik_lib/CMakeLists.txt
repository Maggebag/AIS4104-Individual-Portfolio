find_package(NLopt CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(orocos_kdl CONFIG REQUIRED)

include(FetchContent)
FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG 27cb4c76708608465c413f6d0e6b8d99a4d84302
)
FetchContent_MakeAvailable(spdlog)

set(PUBLIC_HEADERS
    include/trac_ik/math3d.h
    include/trac_ik/kdl_tl.hpp
    include/trac_ik/trac_ik.hpp
    include/trac_ik/nlopt_ik.hpp
    include/trac_ik/dual_quaternion.h
)

add_library(trac_ik
    src/kdl_tl.cpp
    src/nlopt_ik.cpp
    src/trac_ik.cpp
    ${PUBLIC_HEADERS}
)

target_link_libraries(trac_ik PUBLIC NLopt::nlopt spdlog Eigen3::Eigen orocos-kdl)
target_include_directories(trac_ik PUBLIC include/)