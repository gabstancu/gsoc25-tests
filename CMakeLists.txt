cmake_minimum_required(VERSION 3.10)
project(test-project)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include_directories(
    $ENV{HOME}/libs/volesti/include
    $ENV{HOME}/libs/volesti/external
    $ENV{HOMEBREW_OPT_PATH}/eigen/include/eigen3 
    $ENV{HOMEBREW_OPT_PATH}/boost/include
    $ENV{HOMEBREW_OPT_PATH}/lp_solve/include 
    $ENV{HOMEBREW_INCLUDE_PATH} # for mps/mps.h
    $ENV{HOME}/libs/ifopt/ifopt_core/include
    $ENV{HOME}/libs/ifopt/ifopt_ipopt/include
    $ENV{HOMEBREW_OPT_PATH}/gmp/include # for gmp.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include # this project's include statements
)

message("current working directory ${CMAKE_CURRENT_SOURCE_DIR}")

link_directories(
    $ENV{HOMEBREW_OPT_PATH}/boost/lib 
    # /Users/gabrielastancu/libs/lp_solve/lib # lp_solve library path
    $ENV{LIBRARY_PATH}
)

find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_library(LPSOLVE_LIB lpsolve55 PATHS $ENV{HOME}/libs/lp_solve/lib $ENV{LIBRARY_PATH})

add_executable(test-project src/IPM.cpp src/main.cpp)

target_link_libraries(test-project ${LPSOLVE_LIB})