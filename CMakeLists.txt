cmake_minimum_required(VERSION 2.8.12)
project(SCOPE)

include_directories(include/)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_CXX_FLAGS "-std=c++0x -O3")
else()
    set(CMAKE_CXX_FLAGS "-std=c++0x -O3 -lrt -pthread")
endif()

IF (NOT DEFINED SSE_SUPPORT)
    SET(SSE_SUPPORT 0)
ENDIF()
ADD_DEFINITIONS(-DSSE_SUPPORT=${SSE_SUPPORT})

IF (NOT DEFINED DEBUG)
    SET(DEBUG 0)
ENDIF()
ADD_DEFINITIONS(-DDEBUG=${DEBUG})

add_executable(scope src/run_alstructure.cpp src/alstructure.cpp src/matmult.cpp src/genotype.cpp src/storage.cpp)
