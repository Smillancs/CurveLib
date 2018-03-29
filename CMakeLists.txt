cmake_minimum_required(VERSION 3.9)


#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -pthread -Og -ggdb -g3")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -pthread -O3")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")


file(COPY Assets DESTINATION ${CMAKE_BINARY_DIR})

project(Curve_Lib)

add_subdirectory(alglib)
add_subdirectory(CurveLib)
add_subdirectory(GPUcompute)
add_subdirectory(Imgui)
add_subdirectory(Viewer)
#add_subdirectory(Curve_DT)
add_subdirectory(CommandLine)
