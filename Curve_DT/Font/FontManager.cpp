#include "FontManager.h"


FontManager::FontManager()
{
	try
	{
		auto error = FT_Init_FreeType(&library);
		if (error)
		{
			throw std::runtime_error("Error at initializing FreeType library handle.");
		}
	}
	catch (std::exception const &e)
	{
		std::cout << "Exception: " << e.what() << "\n";
	}
}

void FontManager::loadFont(const std::string & ttf, const std::string &name)
{
	faces.insert(std::make_pair(name,FontFace(library, ttf)));
}

FontFace FontManager::getFont(const std::string & ttf)
{
	if (faces.find(ttf)!=faces.end())
		return faces[ttf];
	loadFont(ttf,ttf);
	return getFont(ttf);
}
