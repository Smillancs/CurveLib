add_library(CurveLib STATIC
all_incl.h
all_incl.cpp
bezier.h
bezier.cpp
bezier_gpu.cpp
GLUtils.h
GLUtils.cpp
vertex.vert
fragment.frag
)

include_directories(../glew-2.0.0/include)
include_directories(../OGLPack/include/glm)
find_library(GLEW_LIB libglew32)
find_library(GLU_LIB GlU32)
find_library(OGL_LIB OPENGL32)
target_link_libraries(CurveLib PUBLIC ${GLEW_LIB})
target_link_libraries(CurveLib PUBLIC ${GLU_LIB})
target_link_libraries(CurveLib PUBLIC ${OGL_LIB})
include_directories(../eigen)