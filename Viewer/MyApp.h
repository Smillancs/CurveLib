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
#include "../GPUcompute/gShaderProgram.h"
#include "gVertexBuffer.h"

#include "../CurveLib/Examples.hpp"
#include "../CurveLib/bezierCurve.h"
#include "CurveRenderer.h"

class CMyApp
{
public:
	CMyApp();
	~CMyApp(void);

	bool Init();
  void InitCurveRenderer(Curve::Ptr c, bool gputesselation);
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
	void addCommonShaderAttrib(gShaderProgram&);

	const static int N = 50;
	int vertsInPatch = 2;
protected:
	SDL_Window* win;
	SDL_GLContext context;

	gCamera			m_camera;
	gShaderProgram	m_program_curve;
	gShaderProgram	m_program_tess;
	gShaderProgram	m_program_bez;
	gShaderProgram	m_program_basic;
	gVertexBuffer	m_vb;
	gVertexBuffer axes;

  Curve::Ptr activeCurve;
  int curveSelect = -1;

  Curve::Ptr randomCurve;

  Curve::Ptr optCurve;

	bool tesselated = true;

  bool change = false;

  bool opt = false;
  bool makeOpt = false;
};
