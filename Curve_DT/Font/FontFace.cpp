#include "FontFace.h"

FontFace::FontFace(FT_Library &lib, const std::string &font)
{
	try
	{
		auto error = FT_New_Face(lib, font.c_str(), 0, &face); //todo: change 0 to a proper face index
		if (error)
		{
			throw std::exception("Error at loading a font file!");
		}
	}
	catch (std::exception const &e)
	{
		std::cout << "Exception: " << e.what() << "\n";
	}
}

Glyph FontFace::getGlyph(const wchar_t charCode)
{
	auto glyph_idx = FT_Get_Char_Index(face, charCode);

	auto error = FT_Load_Glyph(face, glyph_idx, FT_LOAD_DEFAULT);

	FT_Glyph glyph;
	error = FT_Get_Glyph(face->glyph, &glyph);

	FT_OutlineGlyph glyphCps = (FT_OutlineGlyph)glyph;

	Glyph bez;

	return bez;
}

void FontFace::setCharSize(const size_t pt)
{
	auto error = FT_Set_Char_Size(
		face,    /* handle to face object           */
		0,       /* char_width in 1/64th of points  */
		pt * 64,   /* char_height in 1/64th of points */
		dpi,     /* horizontal device resolution    */
		dpi);   /* vertical device resolution      */
}
