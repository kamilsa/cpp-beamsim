cmake_minimum_required(VERSION 3.25)

project(ns3-test)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(MPI REQUIRED)

set(ns3_DIR /Users/taisei/dev/beam_ns3_py_0/ns-allinone-3.44/install/lib/cmake/ns3)
find_package(ns3 CONFIG REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main
  ns3::applications
  ns3::core
  ns3::internet
  ns3::point-to-point
)
