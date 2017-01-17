add_executable(Differential_Geometry WIN32
main.cpp
CurveRenderer.h
CurveRenderer.cpp
CommandLine.hpp
gCamera.h
gCamera.cpp
Examples.hpp
Examples.cpp
gShaderProgram.h
gShaderProgram.cpp
gVertexBuffer.h
gVertexBuffer.cpp
MyApp.h
MyApp.cpp
)

set(tarver "x86")


if(UNIX)
	set(libdir "usr/lib/")
	include_directories(/usr/include/SDL2)
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

find_library(SDL2main_LIB NAME SDL2main HINTS ${libdir})
find_library(SDL2_LIB NAME SDL2 HINTS ${libdir})
find_library(SDL2_image_LIB NAME SDL2_image HINTS ${libdir})

target_link_libraries(Differential_Geometry PUBLIC ${GLEW_LIB})
target_link_libraries(Differential_Geometry PUBLIC ${GLU_LIB})
target_link_libraries(Differential_Geometry PUBLIC ${OGL_LIB})
target_link_libraries(Differential_Geometry PUBLIC ${SDL2main_LIB})
target_link_libraries(Differential_Geometry PUBLIC ${SDL2_LIB})
target_link_libraries(Differential_Geometry PUBLIC ${SDL2_image_LIB})
target_link_libraries(Differential_Geometry PUBLIC CurveLib)
