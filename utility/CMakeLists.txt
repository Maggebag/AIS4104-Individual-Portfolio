find_package(Eigen3 CONFIG REQUIRED)
find_package(orocos_kdl CONFIG REQUIRED)

set(PUBLIC_HEADERS
    include/utility/csv.h
    include/utility/enum.h
    include/utility/math.h
    include/utility/vectors.h
)

add_library(utility
    src/math.cpp
    src/vectors.cpp

    ${PUBLIC_HEADERS}
)

target_pedantic_build(utility)
target_link_libraries(utility PUBLIC Eigen3::Eigen trac_ik orocos-kdl)
target_include_directories(utility PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/)