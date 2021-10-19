cmake_minimum_required(VERSION 3.17)
project(optymalizacja)

set(CMAKE_CXX_STANDARD 20)

add_executable(optymalizacja main.cpp matrix.cpp ode_solver.cpp opt_alg.cpp solution.cpp)