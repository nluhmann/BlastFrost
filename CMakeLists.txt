cmake_minimum_required(VERSION 2.8.12)

project(BlastFrost)

#NOTE: If Bifrost is compiled with a different maximum k-mer size, please also adjust this parameter here:
#add_definitions(-DMAX_KMER_SIZE=32)
set(MAX_KMER_SIZE "32")



#set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} "-std=c99 -march=native -Wall -Wextra")
#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++14 -march=native -Wall -Wextra")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -march=native")


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  add_compile_options(-O3)
endif(NOT CMAKE_BUILD_TYPE)


if(CMAKE_BUILD_TYPE MATCHES Debug)
  message("Debug mode")
  add_compile_options(-g)
endif(CMAKE_BUILD_TYPE MATCHES Debug)


#find_library(pthread REQUIRED)
#find_library(z REQUIRED)
#find_library(bifrost REQUIRED)
#find_library(gsl REQUIRED)
#find_library(gslcblas REQUIRED)

file(GLOB SOURCES "src/*.cpp")

#file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp *.hh *.tcc)


# Bifrost library
include(ExternalProject)

ExternalProject_Add(
  bifrost_lib
  BUILD_IN_SOURCE 1
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/bifrost
  PREFIX ${PROJECT_SOURCE_DIR}/bifrost
  CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DMAX_KMER_SIZE=${MAX_KMER_SIZE} -DCOMPILATION_ARCH=OFF DENABLE_AVX2=OFF
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")


ExternalProject_Get_property(bifrost_lib INSTALL_DIR)

add_library(bifrost STATIC IMPORTED)

set(BIFROST_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/bifrost/src")
set(BIFROST_LIBRARY_DIR "${PROJECT_SOURCE_DIR}/bifrost/src")

include_directories(SYSTEM ${BIFROST_INCLUDE_DIR})

set_property(TARGET bifrost PROPERTY IMPORTED_LOCATION "${BIFROST_LIBRARY_DIR}/libbifrost.a")

add_dependencies(bifrost bifrost_lib)

find_library(pthread REQUIRED)
find_library(z REQUIRED)


#file(GLOB SOURCES "src/*.cpp")
add_executable(BlastFrost ${SOURCES})


target_link_libraries(BlastFrost pthread z bifrost)
