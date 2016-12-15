cmake_minimum_required(VERSION 3.6)
 
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")


project(Curve_Lib)

add_subdirectory(CurveLib)
add_subdirectory(Differential_Geometry)
include_directories(eigen)
include_directories(CurveLib)
include_directories(Differential_Geometry)
