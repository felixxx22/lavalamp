cmake_minimum_required(VERSION 3.15)

project(lavalamp)

set(CMAKE_CXX_STANDARD 20)

include(FetchContent)

FetchContent_Declare(
    raylib
    GIT_REPOSITORY https://github.com/raysan5/raylib.git
    GIT_TAG 5.5
)

FetchContent_MakeAvailable(raylib)

add_executable(lavalamp src/stam99.cpp)

target_link_libraries(lavalamp PRIVATE raylib)