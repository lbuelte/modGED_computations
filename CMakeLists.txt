cmake_minimum_required(VERSION 3.3)
project(PolygonClustering)

option(CXX "enable C++ compilation" ON)
if(CXX)
  enable_language(CXX)
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Visual Studio compiler with static runtime libraries
if(MSVC AND MT)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
endif()

#set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Pfad zu Gurobi anpassen
#set(GUROBI_HOME $ENV{GUROBI_HOME})  # Linux Beispiel
#set(GUROBI_HOME "/home/jojo/opt/gurobi1103/linux64")  # Linux Beispiel
#set(GUROBI_HOME "C:/gurobi1000/win64")     # Windows Beispiel

# ignore BOOST deprecated headers
add_definitions("-DBOOST_ALLOW_DEPRECATED_HEADERS")
add_definitions("-DBOOST_BIND_GLOBAL_PLACEHOLDERS")

# end linking pre-built shapelib

find_package(GUROBI REQUIRED)

include_directories(${GUROBI_INCLUDE_DIRS})

include_directories(${GUROBI_HOME}/include)
link_directories(${GUROBI_HOME}/lib)

# # Collect all .cpp files in src/ and all .h files in include/
# file(GLOB SRC_FILES "src/*.cpp")
# file(GLOB INC_FILES "include/*.h")
# set(sources ${SRC_FILES} ${INC_FILES})

# Manually set source files
set(sources
        src/main.cpp
        src/modGED.cpp
        src/graph.cpp
)

add_executable(${CMAKE_PROJECT_NAME} ${sources})

# Inlude paths for private directories
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${GUROBI_INCLUDE_DIRS}
)

# link pre-built shapelib
set(SHAPELIB_DIR "/home/laura/shapelib-1.6.1")  # Adjust if needed
find_library(SHAPELIB_LIB NAMES shp PATHS ${SHAPELIB_DIR}/build)

if (NOT SHAPELIB_LIB)
    message(FATAL_ERROR "shapelib not found")
endif ()

if(CXX)
  set(CMAKE_CXX_STANDARD 23)
  target_link_libraries(${CMAKE_PROJECT_NAME} optimized ${GUROBI_CXX_LIBRARY}
                        debug ${GUROBI_CXX_DEBUG_LIBRARY})
endif()

include_directories(${SHAPELIB_DIR})  # Include shapelib headers
target_link_libraries(${CMAKE_PROJECT_NAME} ${SHAPELIB_LIB})

#target_link_libraries(PolygonClustering gurobi_c++ gurobi110 ${SHAPELIB_LIB})

target_link_libraries(${CMAKE_PROJECT_NAME} ${GUROBI_LIBRARY})

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  include(FeatureSummary)
  feature_summary(WHAT ALL)
endif()
