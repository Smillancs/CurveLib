add_executable(Differential_Geometry WIN32
main.cpp
Curve.hpp
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
GeomInvariant.hpp
GeomInvariant.cpp
MyApp.h
MyApp.cpp
SimpleCurves.hpp
Utils.hpp
)



set(tarver "x86")
set(libdir "../OGLPack/lib/${tarver}")
set(bindir "../OGLPack/bin/${tarver}")


#include_directories(../glew-2.0.0/include)
include_directories(../OGLPack/include/glm)
include_directories(../OGLPack/include)
find_library(GLEW_LIB REQUIRED libglew32 HINTS libdir)
find_library(GLU_LIB REQUIRED GlU32 HINTS libdir)
find_library(OGL_LIB REQUIRED OPENGL32 HINTS libdir)
find_library(SDL2main_LIB NAME SDL2main HINTS libdir)
find_library(SDL2_LIB NAME SDL2 HINTS libdir)
find_library(SDL2_image_LIB NAME SDL2_image HINTS libdir)
target_link_libraries(Differential_Geometry PUBLIC ${GLEW_LIB})
target_link_libraries(Differential_Geometry PUBLIC ${GLU_LIB})
target_link_libraries(Differential_Geometry PUBLIC ${OGL_LIB})
target_link_libraries(Differential_Geometry PUBLIC ${SDL2main_LIB})
target_link_libraries(Differential_Geometry PUBLIC ${SDL2_LIB})
target_link_libraries(Differential_Geometry PUBLIC ${SDL2_image_LIB})
include_directories(../eigen)