add_library(CurveLib STATIC
all_incl.h
all_incl.cpp
bezierCurve.h
bezierCurve.cpp
bezierCurve_gpu.cpp
bezierSquareSurface.h
bezierSquareSurface.cpp
bezierSquareSurface_gpu.cpp
GLUtils.h
GLUtils.cpp
vertex.vert
fragment.frag
Curve.hpp
SimpleCurves.hpp
GeomInvariant.hpp
GeomInvariant.cpp
)


set(tarver "x86")
set(libdir "../OGLPack/lib/${tarver}")
set(bindir "../OGLPack/bin/${tarver}")

include_directories(../eigen)
include_directories(../OGLPack/include)

find_library(GLEW_LIB NAME glew32s HINTS ${libdir})
find_library(OGL_LIB NAME OPENGL32 HINTS ${libdir})
find_library(GLU_LIB NAME GlU32)

target_link_libraries(CurveLib PUBLIC ${GLEW_LIB})
target_link_libraries(CurveLib PUBLIC ${GLU_LIB})
target_link_libraries(CurveLib PUBLIC ${OGL_LIB})