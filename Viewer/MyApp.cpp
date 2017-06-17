#include "MyApp.h"
#include "GLUtils.hpp"

#include "../GPUcompute/GeomOptimize.hpp"
#include "../CurveLib/RandomCurve.hpp"

#include <sstream>

CMyApp::CMyApp()
{
	win = 0;
}

CMyApp::~CMyApp(void)
{
}


bool CMyApp::Init()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);

	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);

  axes = coordAxes();

	std::vector<GeomOptimize::Input2D3> vec = { {glm::vec2(0,0),glm::vec2(1,0),1,-1} };
	GeomOptimize opt;
	std::vector<GeomOptimize::Result> res = opt.optimize2D3("curvatureD", vec);
	BezierCurve optCurve = opt.createResultCurve(vec[0], res[0]);

  Curve::Ptr randomCurve = RandomCurve(5,3);

	//Curve& c = optCurve;
	//Curve& c = ExampleHandler::get(3);
  Curve& c = *randomCurve;

	CurveRenderer ren(c);
	ren.genBufferTesselation(N, 0, 1);
	m_vbT = ren.getBuffer();

	CurveRenderer ren2(c);
	ren2.genBufferNormal(N, 0, 1);
	m_vb = ren2.getBuffer();

	m_program_curve.AttachShader(GL_VERTEX_SHADER, "../Assets/shader_curve.vert");
	m_program_curve.AttachShader(GL_FRAGMENT_SHADER, "../Assets/shader_curve.frag");
	addCommonShaderAttrib(m_program_curve);

	if ( !m_program_curve.LinkProgram() )
	{
		return false;
	}

	m_program_tess.AttachShader(GL_VERTEX_SHADER, "../Assets/shader_tess.vert");
	m_program_tess.AttachShader(GL_TESS_CONTROL_SHADER, "../Assets/beztcs.tcs");
	m_program_tess.AttachShader(GL_TESS_EVALUATION_SHADER, "../Assets/beztes.tes");
	m_program_tess.AttachShader(GL_FRAGMENT_SHADER, "../Assets/shader_tess.frag");
	addCommonShaderAttrib(m_program_tess);

	if ( !m_program_tess.LinkProgram() )
	{
		return false;
	}

	BezierCurve * bez = (BezierCurve *)&ExampleHandler::get(3);
	CurveRenderer ren3(*bez);
	auto cp = bez->GetGlmControlPoints();
	ren3.genBufferCps(cp);
	m_vbBez = ren3.getBuffer();
	vertsInPatch = cp.size();

	m_program_bez.AttachShader(GL_VERTEX_SHADER, "../Assets/shader_tess_adv.vert");
	m_program_bez.AttachShader(GL_TESS_CONTROL_SHADER, "../Assets/shader_tess_adv.tcs");
	m_program_bez.AttachShader(GL_TESS_EVALUATION_SHADER, "../Assets/shader_tess_adv.tes");
	m_program_bez.AttachShader(GL_FRAGMENT_SHADER, "../Assets/shader_tess_adv.frag");
	addCommonShaderAttrib(m_program_bez);

	if ( !m_program_bez.LinkProgram())
	{
		return false;
	}

	m_program_basic.AttachShader(GL_VERTEX_SHADER, "../Assets/shader_basic.vert");
	m_program_basic.AttachShader(GL_FRAGMENT_SHADER, "../Assets/shader_basic.frag");

	m_program_basic.BindAttribLoc(0, "vs_in_pos");
	m_program_basic.BindAttribLoc(1, "vs_in_col");

	if ( !m_program_basic.LinkProgram() )
	{
		return false;
	}

	m_camera.SetProj(45.0f, 640.0f/480.0f, 0.01f, 1000.0f);

	return true;
}


void CMyApp::addCommonShaderAttrib(gShaderProgram& program)
{
	program.BindAttribLoc(0, "vs_in_pos");
	program.BindAttribLoc(1, "vs_in_e");
	program.BindAttribLoc(2, "vs_in_n");
	program.BindAttribLoc(3, "vs_in_b");
	program.BindAttribLoc(4, "vs_in_k");
	program.BindAttribLoc(5, "vs_in_t");
}


void CMyApp::Clean()
{
	m_program_curve.Clean();
}

void CMyApp::Update()
{
	static Uint32 last_time = SDL_GetTicks();
	float delta_time = (SDL_GetTicks() - last_time)/1000.0f;

	m_camera.Update(delta_time);

	last_time = SDL_GetTicks();
}


void CMyApp::Render()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glm::mat4 vp = m_camera.GetViewProj();

	if (tesselated)
	{
		auto & currentProgram = drawBezier ? m_program_bez : m_program_tess;
		currentProgram.On();
		currentProgram.SetUniform("VP", vp);
		if (drawBezier)
		{
			currentProgram.SetUniform( "pointNum", vertsInPatch);
			m_vbBez.On();
			m_vbBez.SetPatchVertices(vertsInPatch);
			m_vbBez.Draw(GL_PATCHES, 0, vertsInPatch);
			m_vbBez.Off();
		}
		else
		{
			m_vbT.On();
			m_vbT.SetPatchVertices(2);
			m_vbT.Draw(GL_PATCHES, 0, 2*(N+1));
			m_vbT.Off();
		}
		currentProgram.Off();
	}
	else
	{
		m_program_curve.On();

		m_program_curve.SetUniform( "VP", vp );

		m_vb.On();
		m_vb.Draw(GL_LINE_STRIP, 0, N+1);
		m_vb.Off();

		m_program_curve.Off();
	}
	m_program_basic.On();

	m_program_basic.SetUniform( "VP", vp );

	axes.On();
	axes.Draw(GL_LINES, 0, 6);
	axes.Off();

	m_program_basic.Off();

}

void CMyApp::KeyboardDown(SDL_KeyboardEvent& key)
{
	m_camera.KeyboardDown(key);
}

void CMyApp::KeyboardUp(SDL_KeyboardEvent& key)
{
	m_camera.KeyboardUp(key);
}

void CMyApp::MouseMove(SDL_MouseMotionEvent& mouse)
{
	m_camera.MouseMove(mouse);
}

void CMyApp::MouseDown(SDL_MouseButtonEvent& mouse)
{
}

void CMyApp::MouseUp(SDL_MouseButtonEvent& mouse)
{
}

void CMyApp::MouseWheel(SDL_MouseWheelEvent& wheel)
{
}

void CMyApp::Resize(int _w, int _h)
{
	glViewport(0, 0, _w, _h);

	m_camera.Resize(_w, _h);
}

void CMyApp::Run()
{
	if (!Init())
	{
		SDL_DestroyWindow(win);
		std::cout << "[app.Init] Az alkalmaz�s inicializ�l�sa k�zben hibat�rt�nt!" << std::endl;
		exit(1);
	}
	bool quit = false;
	SDL_Event ev;
	while (!quit)
	{
		while (SDL_PollEvent(&ev))
		{
			switch (ev.type)
			{
			case SDL_QUIT:
				quit = true;
				break;
			case SDL_KEYDOWN:
				if (ev.key.keysym.sym == SDLK_ESCAPE)
					quit = true;
				KeyboardDown(ev.key);
				break;
			case SDL_KEYUP:
				KeyboardUp(ev.key);
				break;
			case SDL_MOUSEBUTTONDOWN:
				MouseDown(ev.button);
				break;
			case SDL_MOUSEBUTTONUP:
				MouseUp(ev.button);
				break;
			case SDL_MOUSEWHEEL:
				MouseWheel(ev.wheel);
				break;
			case SDL_MOUSEMOTION:
				MouseMove(ev.motion);
				break;
			case SDL_WINDOWEVENT:
				if (ev.window.event == SDL_WINDOWEVENT_SIZE_CHANGED)
				{
					Resize(ev.window.data1, ev.window.data2);
				}
				break;
			}
		}

		Update();
		Render();

		SDL_GL_SwapWindow(win);
	}

	Clean();

	SDL_GL_DeleteContext(context);
	SDL_DestroyWindow(win);

}


void CMyApp::initGraphics()
{
	if (SDL_Init(SDL_INIT_VIDEO) == -1)
	{
		std::cout << "[SDL ind�t�sa]Hiba az SDL inicializ�l�sa k�zben: " << SDL_GetError() << std::endl;
		exit(1);
	}

	SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

	SDL_GL_SetAttribute(SDL_GL_BUFFER_SIZE, 32);
	SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
	SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);

	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

	win = SDL_CreateWindow("Hello SDL&OpenGL!",	100, 100, 640, 480,	SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);

	if (win == 0)
	{
		std::cout << "[Ablak l�trehoz�sa]Hiba az SDL inicializ�l�sa k�zben: " << SDL_GetError() << std::endl;
		exit(1);
	}

	context = SDL_GL_CreateContext(win);
	if (context == 0)
	{
		std::cout << "[OGL context l�trehoz�sa]Hiba az SDL inicializ�l�sa k�zben: " << SDL_GetError() << std::endl;
		exit(1);
	}

	SDL_GL_SetSwapInterval(1);

	GLenum error = glewInit();
	if (error != GLEW_OK)
	{
		std::cout << "[GLEW] Hiba az inicializ�l�s sor�n!" << std::endl;
		exit(1);
	}

	int glVersion[2] = { -1, -1 };
	glGetIntegerv(GL_MAJOR_VERSION, &glVersion[0]);
	glGetIntegerv(GL_MINOR_VERSION, &glVersion[1]);
	std::cout << "Running OpenGL " << glVersion[0] << "." << glVersion[1] << std::endl;

	if (glVersion[0] == -1 && glVersion[1] == -1)
	{
		SDL_GL_DeleteContext(context);
		SDL_DestroyWindow(win);

		std::cout << "[OGL context l�trehoz�sa] Nem siker�lt l�trehozni az OpenGL context-et! Lehet, hogy az SDL_GL_SetAttribute(...) h�v�sokn�l az egyik be�ll�t�s helytelen." << std::endl;

		exit(1);
	}

	std::stringstream window_title;
	window_title << "OpenGL " << glVersion[0] << "." << glVersion[1];
	SDL_SetWindowTitle(win, window_title.str().c_str());

}
