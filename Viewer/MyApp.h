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

class CMyApp
{
public:
	CMyApp();
	~CMyApp(void);

	bool Init();
	void Clean();

	void Update();
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

	const static int N = 50;
protected:
	SDL_Window* win;
	SDL_GLContext context;

	gCamera			m_camera;
	gShaderProgram	m_program_curve;
	gShaderProgram	m_program_tess;
	gShaderProgram	m_program_basic;
	gVertexBuffer	m_vb;
	gVertexBuffer	m_vbT;
	gVertexBuffer   axes;
	
	bool tesselated = true;
};

