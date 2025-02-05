cmake_minimum_required(VERSION 3.5.0)
project(sign-approximation VERSION 0.1.0 LANGUAGES C CXX)

find_library(NTL_LIB NAMES ntl PATHS /usr/lib /usr/local/lib REQUIRED)
find_path(NTL_INCLUDE_DIR NAMES NTL/RR.h PATHS /usr/include /usr/local/include REQUIRED)

include_directories(${NTL_INCLUDE_DIR})

add_executable(sign-approximation remez.cpp Point.cpp)

target_link_libraries(sign-approximation ${NTL_LIB})
