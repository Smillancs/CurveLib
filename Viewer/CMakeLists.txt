add_executable(Viewer WIN32
main.cpp
CurveRenderer.h
CurveRenderer.cpp
gCamera.h
gCamera.cpp
gVertexBuffer.h
gVertexBuffer.cpp
MyApp.h
MyApp.cpp
imgui_impl_sdl_gl3.h
imgui_impl_sdl_gl3.cpp
)


set(tarver "x86")


if(UNIX)
	set(libdir "usr/lib/")
	include_directories(/usr/include/SDL2)
	include_directories(../Imgui)
	find_library(GLEW_LIB REQUIRED GLEW HINTS ${libdir})
	find_library(GLU_LIB REQUIRED GLU HINTS ${libdir})
	find_library(OGL_LIB REQUIRED GL HINTS ${libdir})
endif()
if(WIN32)
	set(libdir "../OGLPack/lib/${tarver}")
	include_directories(../eigen)
	include_directories(../OGLPack/include)
	include_directories(../imgui)
	find_library(GLEW_LIB NAME glew32s HINTS ${libdir})
	find_library(GLU_LIB NAME GlU32 HINTS ${libdir})
	find_library(OGL_LIB NAME OPENGL32 HINTS ${libdir})
endif()

find_library(SDL2main_LIB NAME SDL2main HINTS ${libdir})
find_library(SDL2_LIB NAME SDL2 HINTS ${libdir})
find_library(SDL2_image_LIB NAME SDL2_image HINTS ${libdir})

target_link_libraries(Viewer PUBLIC ${GLEW_LIB})
target_link_libraries(Viewer PUBLIC ${GLU_LIB})
target_link_libraries(Viewer PUBLIC ${OGL_LIB})
target_link_libraries(Viewer PUBLIC ${SDL2main_LIB})
target_link_libraries(Viewer PUBLIC ${SDL2_LIB})
target_link_libraries(Viewer PUBLIC ${SDL2_image_LIB})
target_link_libraries(Viewer PUBLIC CurveLib)
target_link_libraries(Viewer PUBLIC GPUcompute)
target_link_libraries(Viewer PUBLIC Imgui)
