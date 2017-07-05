#pragma once

#include <stdexcept>
#include <iostream>
#include <map>
#include "FontFace.h"

#include <ft2build.h>
#include FT_FREETYPE_H

class FontManager
{
	FT_Library  library;
	std::map<std::string,FontFace> faces;
public:
	FontManager();
	void loadFont(const std::string &ttf, const std::string &name);
	FontFace getFont(const std::string &ttf);
};