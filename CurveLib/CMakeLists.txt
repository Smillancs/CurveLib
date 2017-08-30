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
Curve.hpp
SimpleCurves.hpp
GeomInvariant.hpp
GeomInvariant.cpp
Utils.hpp
Examples.hpp
Examples.cpp
RandomCurve.hpp
RandomCurve.cpp
MathFunctions.hpp
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
	set(bindir "../OGLPack/bin/${tarver}")
	include_directories(../eigen)
	include_directories(../OGLPack/include)

	find_library(GLEW_LIB NAME glew32s HINTS ${libdir})
	find_library(OGL_LIB NAME OPENGL32 HINTS ${libdir})
	find_library(GLU_LIB NAME GlU32)
endif()


target_link_libraries(CurveLib PUBLIC ${GLEW_LIB})
target_link_libraries(CurveLib PUBLIC ${GLU_LIB})
target_link_libraries(CurveLib PUBLIC ${OGL_LIB})
