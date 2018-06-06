cmake_minimum_required(VERSION 3.0)

project (alamo)

# This prevents the user from building inside the source directory
# (which would overwrite the legacy makefile)
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed. You must create another build directory.\nIMPORTANT: be sure to remove CMakeCache.txt and CMakeFiles/ with \nrm -rf CMakeCache.txt CMakeFiles")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# C++11 is required for compatibility with AMREX
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "-Wpedantic -Wextra -Wall -std=c++11 -ggdb")

# Locate source files
file(GLOB_RECURSE alamo_SOURCES "src/*.cpp")
file(GLOB_RECURSE alamo_HEADERS "src/*.h")

# Locate AMReX
find_library(AMREX NAMES "libamrex.a" NO_DEFAULT_PATH)
if(NOT AMREX)
  message(FATAL_ERROR "AMReX library not found\n(Specify location using -DAMREX=/path/to/amrex or running ccmake)")
endif(NOT AMREX)
# Locate MPI
find_package(MPI REQUIRED)

set (SPACEDIM 3 CACHE STRING "Dimensionality of the system")

# Include and Link
include_directories(SYSTEM ${AMREX}/include/ ${MPI_INCLUDE_PATH})
include_directories("src/")
link_directories(${AMREX}/lib/)
link_libraries(${MPI_LIBRARIES} amrex gfortran)


# Determine git hash, username, etc for metadata
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE METADATA_GITHASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
  COMMAND whoami
  OUTPUT_VARIABLE METADATA_USER
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
  COMMAND hostname
  OUTPUT_VARIABLE METADATA_HOSTNAME
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
  COMMAND date +%x
  OUTPUT_VARIABLE METADATA_DATE
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
  COMMAND date +%H:%M:%S
  OUTPUT_VARIABLE METADATA_TIME
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_definitions("-DMETADATA_GITHASH=\"${METADATA_GITHASH}\"")
add_definitions("-DMETADATA_USER=\"${METADATA_USER}\"")
add_definitions("-DMETADATA_PLATFORM=\"${METADATA_PLATFORM}\"")
add_definitions("-DMETADATA_COMPILER=\"${METADATA_COMPILER}\"")
add_definitions("-DMETADATA_DATE=\"${METADATA_DATE}\"")
add_definitions("-DMETADATA_TIME=\"${METADATA_TIME}\"")
add_definitions("-DAMREX_SPACEDIM=${SPACEDIM}")

# Compile source
add_library(alamo_lib OBJECT ${alamo_SOURCES})

# Specify Executables
add_executable (alamo "src/alamo.cc" $<TARGET_OBJECTS:alamo_lib>)
add_executable (fem "src/fem.cc" $<TARGET_OBJECTS:alamo_lib>)
add_executable (flame "src/flame.cc" $<TARGET_OBJECTS:alamo_lib>)
add_executable (degradation "src/degradation.cc" $<TARGET_OBJECTS:alamo_lib>)

