
function(add_unit_test filename)
    get_filename_component(name ${filename} NAME_WE)
    add_executable(${name} ${name}.cpp)
    target_pedantic_build(${name})
    target_link_libraries(${name} PRIVATE Catch2::Catch2WithMain)
    target_include_directories(${name} PRIVATE "${PROJECT_SOURCE_DIR}/src")
    add_test(NAME ${name} COMMAND ${name})
endfunction()

add_subdirectory(utility)
