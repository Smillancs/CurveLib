cmake_minimum_required(VERSION 3.5)
 
set(CMAKE_CXX_FLAGS "-std=c++11 -g")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")


project(Curve_Lib)


add_subdirectory(CurveLib)
add_subdirectory(GPUcompute)
add_subdirectory(Imgui)
add_subdirectory(Viewer)
add_subdirectory(CommandLine)
