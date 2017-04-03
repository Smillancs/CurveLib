cmake_minimum_required(VERSION 3.5)
 
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")


project(Curve_Lib)

add_subdirectory(CurveLib)
add_subdirectory(Imgui)
add_subdirectory(Viewer)
add_subdirectory(Curve_DT)
add_subdirectory(CommandLine)
