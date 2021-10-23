cmake_minimum_required(VERSION 3.17)
project(metody-optymalizacji)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g -ggdb -Wall -trigraphs")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-Wall")

set(SOURCES main.cpp matrix.cpp ode_solver.cpp opt_alg.cpp solution.cpp)

add_executable(metody-optymalizacji ${SOURCES})