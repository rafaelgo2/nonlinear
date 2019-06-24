cmake_minimum_required(VERSION 3.5)
project(tp2)

SET(CMAKE_C_FLAGS_DEBUG "-D_DEBUG")

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3")

#find_package(Armadillo REQUIRED)
#include_directories(${ARMADILLO_INCLUDE_DIRS})

set(SOURCE_FILES main.cpp)
add_executable(tp2 ${SOURCE_FILES})
#target_link_libraries(spec ${ARMADILLO_LIBRARIES})
