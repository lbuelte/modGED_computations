cmake_minimum_required(VERSION 3.30)
project(PolygonClustering)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Pfad zu Gurobi anpassen
#set(GUROBI_HOME $ENV{GUROBI_HOME})  # Linux Beispiel
set(GUROBI_HOME "/home/jojo/opt/gurobi1103/linux64")  # Linux Beispiel
#set(GUROBI_HOME "C:/gurobi1000/win64")     # Windows Beispiel

# ignore BOOST deprecated headers
add_definitions("-DBOOST_ALLOW_DEPRECATED_HEADERS")
add_definitions("-DBOOST_BIND_GLOBAL_PLACEHOLDERS")

# end linking pre-built shapelib

include_directories(${GUROBI_HOME}/include)
link_directories(${GUROBI_HOME}/lib)

# source files
set(sources
        src/main.cpp
        src/distance_lin.cpp
        include/distance_lin.h
        include/binary_io.h
        include/graph_boost.h
        src/binary_io.cpp
        src/graph.cpp
        include/shapefile_io_operations.h
        src/shapefile_io_operations.cpp
        include/hierachicalClustering3.h
        src/hierachicalClustering3.cpp
)

add_executable(${CMAKE_PROJECT_NAME} ${sources})

# link pre-built shapelib
set(SHAPELIB_DIR "${CMAKE_SOURCE_DIR}/extern/shapelib")  # Adjust if needed
find_library(SHAPELIB_LIB NAMES shp PATHS ${SHAPELIB_DIR}/build)

if (NOT SHAPELIB_LIB)
    message(FATAL_ERROR "shapelib not found")
endif ()

include_directories(${SHAPELIB_DIR})  # Include shapelib headers
target_link_libraries(${CMAKE_PROJECT_NAME} ${SHAPELIB_LIB})

target_link_libraries(PolygonClustering gurobi_c++ gurobi110 ${SHAPELIB_LIB})
