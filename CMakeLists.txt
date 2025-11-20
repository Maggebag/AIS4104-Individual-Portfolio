cmake_minimum_required(VERSION 3.30)
project(AIS4104-Portfolio)

set(CMAKE_CXX_STANDARD 20)

option(BUILD_TESTS "Build unit tests" ON)
option(BUILD_ALL_WARNINGS "Build with all compiler warnings" OFF)

function(target_pedantic_build TARGET)
    if (${BUILD_ALL_WARNINGS})
        target_compile_options(${TARGET} PRIVATE
            $<$<CXX_COMPILER_ID:MSVC>:/W4>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>
        )
        target_compile_definitions(${TARGET} PRIVATE _CRT_SECURE_NO_WARNINGS)
    endif ()
endfunction()

add_subdirectory(third_party)

add_subdirectory(app)
add_subdirectory(utility)
add_subdirectory(simulation)
add_subdirectory(visualization)

add_dependencies(app utility)
add_dependencies(app simulation)
add_dependencies(app visualization)
add_dependencies(simulation utility)
add_dependencies(visualization utility)

# Deploy the URDF and models folders to the app runtime directory
add_custom_target(install_model_files)
add_dependencies(install_model_files app)
add_custom_command(TARGET install_model_files PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/third_party/urdf/ $<TARGET_FILE_DIR:app>/urdf)
add_custom_command(TARGET install_model_files PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/third_party/models/ $<TARGET_FILE_DIR:app>/models)

if (BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.4.0
    )
    FetchContent_MakeAvailable(Catch2)

    enable_testing()
    add_subdirectory(tests)
endif ()