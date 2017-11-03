add_executable(CommandLine WIN32
main.cpp
CommandLine.hpp
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
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /subsystem:console /ENTRY:mainCRTStartup")
	set(libdir "../OGLPack/lib/${tarver}")
	include_directories(../eigen)
	include_directories(../OGLPack/include)
	find_library(GLEW_LIB NAME glew32 HINTS ${libdir})
	find_library(GLU_LIB NAME GlU32 HINTS ${libdir})
	find_library(OGL_LIB NAME opengl32)
endif()

find_library(SDL2main_LIB NAME SDL2main HINTS ${libdir})
find_library(SDL2_LIB NAME SDL2 HINTS ${libdir})
find_library(SDL2_image_LIB NAME SDL2_image HINTS ${libdir})

target_link_libraries(CommandLine PUBLIC ${GLEW_LIB})
if(UNIX)
	target_link_libraries(CommandLine PUBLIC ${GLU_LIB})
	target_link_libraries(CommandLine PUBLIC ${OGL_LIB})
endif()
if(WIN32)
	SET_TARGET_PROPERTIES(CommandLine PROPERTIES LINK_FLAGS "/link opengl32.lib")
endif()
target_link_libraries(CommandLine PUBLIC ${SDL2_LIB})
target_link_libraries(CommandLine PUBLIC CurveLib)
target_link_libraries(CommandLine PUBLIC GPUcompute)


if(WIN32)
	file(COPY ../OGLPack/bin/x86/glew32.dll DESTINATION ${CMAKE_BINARY_DIR}/CommandLine)
	file(COPY ../OGLPack/bin/x86/SDL2.dll DESTINATION ${CMAKE_BINARY_DIR}/CommandLine)
endif()