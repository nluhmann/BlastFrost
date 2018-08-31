cmake_minimum_required(VERSION 2.8.12)

project(BlastFrost)


set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-std=c99 -march=native")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++14 -march=native")


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  add_compile_options(-O3)
endif(NOT CMAKE_BUILD_TYPE)


if(CMAKE_BUILD_TYPE MATCHES Debug)
  message("Debug mode")
  add_compile_options(-g)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

find_library(pthread REQUIRED)
find_library(z REQUIRED)
find_library(bifrost REQUIRED)
find_library(gsl REQUIRED)
find_library(gslcblas REQUIRED)

file(GLOB SOURCES "src/*.cpp")
add_executable(BlastFrost ${SOURCES})


target_link_libraries(BlastFrost pthread z bifrost gsl gslcblas)
