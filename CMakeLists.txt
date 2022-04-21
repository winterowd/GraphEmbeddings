# polyakov loop effective theories project
project(PLEFT)

cmake_minimum_required(VERSION 2.8.12)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

###### options ######

option(USE_DEBUG "Enter debug mode" OFF)
if (USE_DEBUG)
    add_definitions( -DDEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

###### External Libraries ######

# Boost
find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Boost_LIBRARIES})

# CLN
find_package(CLN COMPONENTS program_options REQUIRED)
include_directories(${CLN_INCLUDE_DIRS})
set(LIBS ${LIBS} ${CLN_LIBRARIES})

# GiNaC (by hand from user)
if (DEFINED GINAC_DIR)
    message(STATUS "GiNaC located at ${GINAC_DIR}!")
else()
    set(GINAC_DIR "/usr/local")
    message(STATUS "GiNaC directory not supplied by user! (default used: ${GINAC_DIR})")
endif()
add_library(ginac STATIC IMPORTED)
set_target_properties(ginac PROPERTIES IMPORTED_LOCATION ${GINAC_DIR}/lib/libginac.a)
include_directories(${GINAC_DIR}/include/ginac/)

# Nauty (by hand from user)
if (DEFINED NAUTY_DIR)
   message(STATUS "Nauty located at ${NAUTY_DIR}!")
else()
   message(FATAL_ERROR "Nauty directory not supplied by user!")
endif()
add_library(nauty STATIC IMPORTED)
set_target_properties(nauty PROPERTIES IMPORTED_LOCATION ${NAUTY_DIR}/nauty.a)
include_directories(${NAUTY_DIR})

###### Internal Libraries ######

#include Graph
add_subdirectory(Graph)
include_directories(Graph)

set(PROJECT_LIBS ${PROJECT_LIBS} Graph)

###### all executables ######

# Embed
add_executable(Embed Programs/Embed.cpp)
target_link_libraries(Embed ${PROJECT_LIBS} ${LIBS} nauty ginac)

# GenerateNauty
add_executable(GenerateNauty ${NAUTY_OBJS} Programs/GenerateNauty.cpp)
target_link_libraries(GenerateNauty ${PROJECT_LIBS} ${LIBS} nauty)

###### all tests ######

add_executable(TestCubicLatticeCanonicalizor Tests/TestCubicLatticeCanonicalizor.cpp)
target_link_libraries(TestCubicLatticeCanonicalizor ${PROJECT_LIBS} ${LIBS} nauty ginac)

add_executable(TestGraphWeight Tests/TestGraphWeight.cpp)
target_link_libraries(TestGraphWeight ${PROJECT_LIBS} ${LIBS} nauty ginac)

add_executable(TestZCluster Tests/TestZCluster.cpp)
target_link_libraries(TestZCluster ${PROJECT_LIBS} ${LIBS} nauty ginac)

add_executable(TestPowerSet Tests/TestPowerSet.cpp)
target_link_libraries(TestPowerSet ${PROJECT_LIBS} nauty)

add_executable(TestDenseNauty Tests/TestDenseNauty.cpp)
target_link_libraries(TestDenseNauty ${PROJECT_LIBS} nauty)

add_executable(TestRecursionXi Tests/TestRecursionXi.cpp)
target_link_libraries(TestRecursionXi ${PROJECT_LIBS} ${LIBS} nauty ginac)

add_executable(TestGinacPolynomial Tests/TestGinacPolynomial.cpp)
target_link_libraries(TestGinacPolynomial ${PROJECT_LIBS} ${LIBS} nauty ginac)

add_executable(TestCorrelator Tests/TestCorrelator.cpp)
target_link_libraries(TestCorrelator ${PROJECT_LIBS} ${LIBS} nauty ginac)
