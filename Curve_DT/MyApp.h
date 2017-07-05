#pragma once

// GLEW
#include <GL/glew.h>

// SDL
#include <SDL.h>
#include <SDL_opengl.h>

// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform2.hpp>
#include <glm/gtc/constants.hpp>

#include "gCamera.h"
#include "gShaderProgram.h"
#include "gVertexBuffer.h"

#include "../CurveLib/Examples.hpp"
#include "../CurveLib/bezierCurve.h"
#include "CurveRenderer.h"

#include "Font/FontManager.h"
#include "Font/Glyph.h"

class CMyApp
{
public:
	CMyApp();
	~CMyApp(void);

	bool Init();
	void Clean();

	void Update();
	glm::vec3 QuerryAndShow(const std::vector<glm::vec3> & cp, const int size, const glm::vec3 &point, const glm::mat4 &trf);
	void DrawHull(const gVertexBuffer & from, const int size, const glm::mat4 &trf);
	void DrawGuess(const glm::vec3 &point, const glm::vec3 &guess);
	void Render();

	void KeyboardDown(SDL_KeyboardEvent&);
	void KeyboardUp(SDL_KeyboardEvent&);
	void MouseMove(SDL_MouseMotionEvent&);
	void MouseDown(SDL_MouseButtonEvent&);
	void MouseUp(SDL_MouseButtonEvent&);
	void MouseWheel(SDL_MouseWheelEvent&);
	void Resize(int, int);

	void Run();

	void initGraphics();

private:
	const static int N = 50;
	int vertsInPatch = 2;

	std::vector<glm::vec3> glmControlPoints;
	Curve* currentCurve = &ExampleHandler::get(3);
	CurveRenderer* generator;

	GLuint fb, colorBuffer, pointBuffer;
	
	FontManager f;
	void SetUpDummyFont();

protected:
	SDL_Window* win;
	SDL_GLContext context;
	
	glm::vec2 res;
	gCamera			m_camera;
	gShaderProgram	m_program_bez, m_program_lines;
	gVertexBuffer	m_vbBez, m_vbQuad;	
};

