cmake_minimum_required(VERSION 3.30)
project(ngmToRoot)

find_package(ROOT REQUIRED)

add_definitions(${ROOT_DEFINITIONS})
set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
include_directories(${ROOT_INCLUDE_DIRS})
link_libraries(${ROOT_LIBRARIES})

add_executable(ngmToRoot main.cpp)
