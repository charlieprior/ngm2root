cmake_minimum_required(VERSION 3.30)
project(ngm2root)

add_executable(${PROJECT_NAME} main.cpp)

# fetch latest argparse
include(FetchContent)
FetchContent_Declare(
        argparse
        GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

## for boost
find_package(Boost 1.67 REQUIRED) # sort is header only. 1.67 is minimum version for spinsort I think
target_include_directories(${PROJECT_NAME} PRIVATE ${Boost_INCLUDE_DIRS})

# for ROOT
find_package(ROOT REQUIRED)
add_definitions(${ROOT_DEFINITIONS}) # unfortunately required by root
set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD}) # unfortunately required by root
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}") # unfortunately required by root
target_include_directories(${PROJECT_NAME} PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${ROOT_LIBRARIES})

target_precompile_headers(${PROJECT_NAME} PRIVATE <boost/sort/sort.hpp> <argparse/argparse.hpp>)
target_link_libraries(${PROJECT_NAME} PRIVATE argparse)

target_compile_options(${PROJECT_NAME} PRIVATE
        $<$<CXX_COMPILER_ID:MSVC>:/W4>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>
)