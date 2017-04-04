add_library(GPUcompute STATIC
gShaderProgram.h
gShaderProgram.cpp
gBuffer.h
gBuffer.cpp
GeomOptimize.hpp
GeomOptimize.cpp
)


set(tarver "x86")


if(UNIX)
	set(libdir "usr/lib/")
	find_library(GLEW_LIB REQUIRED GLEW HINTS ${libdir})
	find_library(GLU_LIB REQUIRED GLU HINTS ${libdir})
	find_library(OGL_LIB REQUIRED GL HINTS ${libdir})
endif()
if(WIN32)
	set(libdir "../OGLPack/lib/${tarver}")
	include_directories(../eigen)
	include_directories(../OGLPack/include)
	find_library(GLEW_LIB NAME glew32s HINTS ${libdir})
	find_library(GLU_LIB NAME GlU32 HINTS ${libdir})
	find_library(OGL_LIB NAME OPENGL32 HINTS ${libdir})
endif()

target_link_libraries(GPUcompute PUBLIC ${GLEW_LIB})
target_link_libraries(GPUcompute PUBLIC ${GLU_LIB})
target_link_libraries(GPUcompute PUBLIC ${OGL_LIB})
target_link_libraries(GPUcompute PUBLIC CurveLib)
