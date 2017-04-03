add_executable(Curve_DT WIN32
main.cpp
CurveRenderer.h
CurveRenderer.cpp
gCamera.h
gCamera.cpp
gShaderProgram.h
gShaderProgram.cpp
gVertexBuffer.h
gVertexBuffer.cpp
MyApp.h
MyApp.cpp
imgui_impl_sdl_gl3.h
imgui_impl_sdl_gl3.cpp
${CMAKE_BINARY_DIR}/Assets/BezierEval.glsl
${CMAKE_BINARY_DIR}/Assets/GeomInvariant.glsl
${CMAKE_BINARY_DIR}/Assets/shader_basic.frag
${CMAKE_BINARY_DIR}/Assets/shader_basic.vert
${CMAKE_BINARY_DIR}/Assets/shader_tess_adv.frag
${CMAKE_BINARY_DIR}/Assets/shader_tess_adv.tcs
${CMAKE_BINARY_DIR}/Assets/shader_tess_adv.tes
${CMAKE_BINARY_DIR}/Assets/shader_tess_adv.vert
)

source_group( Shaders FILES
${CMAKE_BINARY_DIR}/Assets/BezierEval.glsl
${CMAKE_BINARY_DIR}/Assets/GeomInvariant.glsl
${CMAKE_BINARY_DIR}/Assets/shader_basic.frag
${CMAKE_BINARY_DIR}/Assets/shader_basic.vert
${CMAKE_BINARY_DIR}/Assets/shader_tess_adv.frag
${CMAKE_BINARY_DIR}/Assets/shader_tess_adv.tcs
${CMAKE_BINARY_DIR}/Assets/shader_tess_adv.tes
${CMAKE_BINARY_DIR}/Assets/shader_tess_adv.vert
)

set(tarver "x86")

if(UNIX)
	set(libdir "usr/lib/")
	include_directories(/usr/include/SDL2)
	include_directories(../Imgui)
	include_directories(../Assets)
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

target_link_libraries(Curve_DT PUBLIC ${GLEW_LIB})
target_link_libraries(Curve_DT PUBLIC ${GLU_LIB})
target_link_libraries(Curve_DT PUBLIC ${OGL_LIB})
target_link_libraries(Curve_DT PUBLIC ${SDL2main_LIB})
target_link_libraries(Curve_DT PUBLIC ${SDL2_LIB})
target_link_libraries(Curve_DT PUBLIC ${SDL2_image_LIB})
target_link_libraries(Curve_DT PUBLIC CurveLib)
target_link_libraries(Curve_DT PUBLIC Imgui)