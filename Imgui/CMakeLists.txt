add_library(Imgui STATIC
imconfig.h
imgui.cpp
imgui.h
imgui_demo.cpp
imgui_draw.cpp
imgui_internal.h
stb_rect_pack.h
stb_textedit.h
stb_truetype.h
)

set(tarver "x86")


if(UNIX)
	set(libdir "usr/lib/")
endif()
if(WIN32)
	set(libdir "../OGLPack/lib/${tarver}")
endif()