cmake_minimum_required(VERSION 3.16)

project(PressureVelocityCoupling LANGUAGES CXX)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# including the folders as libraries
add_library(mesh INTERFACE)
target_include_directories(mesh INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Mesh
    )

add_library(sparse_operations INTERFACE)
target_include_directories(sparse_operations INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/SparseOperations
    )
target_link_libraries(sparse_operations INTERFACE mesh)

add_library(solver INTERFACE)
target_include_directories(solver INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Solver
    )
target_link_libraries(solver INTERFACE sparse_operations)

add_library(fvcore INTERFACE)
target_include_directories(fvcore INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/FVCore
)
target_link_libraries(fvcore INTERFACE mesh sparse_operations)

add_library(discretisation INTERFACE)
target_include_directories(discretisation INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Discretisation
)
target_link_libraries(discretisation INTERFACE fvcore)

add_library(test INTERFACE)
target_include_directories(test INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Test
)
target_link_libraries(test INTERFACE fvcore mesh sparse_operations)

add_executable(pvc_app main.cpp)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(pvc_app PRIVATE
        -Wall
        -Wextra
        -Wpedantic
    )
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(pvc_app PRIVATE
        /W4
    )
endif()

target_link_libraries(pvc_app PRIVATE
    mesh
    sparse_operations
    fvcore
    discretisation
    solver
    test
)

add_custom_command(TARGET pvc_app POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_CURRENT_SOURCE_DIR}/Files
            $<TARGET_FILE_DIR:pvc_app>/Files
)

add_custom_command(TARGET pvc_app POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
            $<TARGET_FILE_DIR:pvc_app>/Output
)
