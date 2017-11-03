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
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /subsystem:console /ENTRY:mainCRTStartup")
	set(libdir "../OGLPack/lib/${tarver}")
	include_directories(../eigen)
	include_directories(../OGLPack/include)
	include_directories(../imgui)
	find_library(GLEW_LIB NAME glew32 HINTS ${libdir})
	find_library(GLU_LIB NAME GlU32 HINTS ${libdir})
	find_library(OGL_LIB NAME opengl32)
endif()

find_library(SDL2main_LIB NAME SDL2main HINTS ${libdir})
find_library(SDL2_LIB NAME SDL2 HINTS ${libdir})
find_library(SDL2_image_LIB NAME SDL2_image HINTS ${libdir})

target_link_libraries(Viewer PUBLIC ${GLEW_LIB})
if(UNIX)
	target_link_libraries(Viewer PUBLIC ${GLU_LIB})
	target_link_libraries(Viewer PUBLIC ${OGL_LIB})
endif()
if(WIN32)
	SET_TARGET_PROPERTIES(Viewer PROPERTIES LINK_FLAGS "/link opengl32.lib")
endif()
target_link_libraries(Viewer PUBLIC ${SDL2main_LIB})
target_link_libraries(Viewer PUBLIC ${SDL2_LIB})
target_link_libraries(Viewer PUBLIC ${SDL2_image_LIB})
target_link_libraries(Viewer PUBLIC CurveLib)
target_link_libraries(Viewer PUBLIC GPUcompute)
target_link_libraries(Viewer PUBLIC Imgui)
SET_TARGET_PROPERTIES(Viewer PROPERTIES LINK_FLAGS "/link opengl32.lib")

if(WIN32)
	file(COPY ../OGLPack/bin/x86/glew32.dll DESTINATION ${CMAKE_BINARY_DIR}/Viewer)
	file(COPY ../OGLPack/bin/x86/SDL2.dll DESTINATION ${CMAKE_BINARY_DIR}/Viewer)
endif()