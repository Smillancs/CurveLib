#pragma once

#include <iostream>
#include <string>

#include "Glyph.h"

#include <ft2build.h>
#include "freetype\ftglyph.h"

class FontFace
{
	FT_Face face;
	const size_t dpi = 300;
public:
	FontFace() {}
	FontFace(FT_Library &lib, const std::string &font);
	
	void setCharSize(const size_t pt);

	Glyph getGlyph(const wchar_t charCode);
};