#include "MyApp.h"
#include "GLUtils.hpp"

#include "../CurveLib/Examples.hpp"
#include "CurveRenderer.h"

//#include <GL/GLU.h>
//#include <math.h>
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

	Curve& c = ExampleHandler::get(3);

	CurveRenderer ren(c);
	ren.genBufferTesselation(N, 0, 1);
	m_vbT = ren.getBuffer();

	CurveRenderer ren2(c);
	ren2.genBufferNormal(N, 0, 1);
	m_vb = ren2.getBuffer();

	m_program_curve.AttachShader(GL_VERTEX_SHADER, "shader_curve.vert");
	m_program_curve.AttachShader(GL_FRAGMENT_SHADER, "shader_curve.frag");

	m_program_curve.BindAttribLoc(0, "vs_in_pos");
	m_program_curve.BindAttribLoc(1, "vs_in_e");
	m_program_curve.BindAttribLoc(2, "vs_in_n");
	m_program_curve.BindAttribLoc(3, "vs_in_b");
	m_program_curve.BindAttribLoc(4, "vs_in_k");
	m_program_curve.BindAttribLoc(5, "vs_in_t");

	if ( !m_program_curve.LinkProgram() )
	{
		return false;
	}

	m_program_tess.AttachShader(GL_VERTEX_SHADER, "shader_tess.vert");
	m_program_tess.AttachShader(GL_TESS_CONTROL_SHADER, "beztcs.tcs");
	m_program_tess.AttachShader(GL_TESS_EVALUATION_SHADER, "beztes.tes");
	m_program_tess.AttachShader(GL_FRAGMENT_SHADER, "shader_tess.frag");

	m_program_tess.BindAttribLoc(0, "vs_in_pos");
	m_program_tess.BindAttribLoc(1, "vs_in_e");
	m_program_tess.BindAttribLoc(2, "vs_in_n");
	m_program_tess.BindAttribLoc(3, "vs_in_b");
	m_program_tess.BindAttribLoc(4, "vs_in_k");
	m_program_tess.BindAttribLoc(5, "vs_in_t");

	if ( !m_program_tess.LinkProgram() )
	{
		return false;
	}

	m_program_basic.AttachShader(GL_VERTEX_SHADER, "shader_basic.vert");
	m_program_basic.AttachShader(GL_FRAGMENT_SHADER, "shader_basic.frag");

	m_program_basic.BindAttribLoc(0, "vs_in_pos");
	m_program_basic.BindAttribLoc(1, "vs_in_col");

	if ( !m_program_basic.LinkProgram() )
	{
		return false;
	}

	m_camera.SetProj(45.0f, 640.0f/480.0f, 0.01f, 1000.0f);

	return true;
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
	if (tesselated)
	{
		m_program_tess.On();

		glm::mat4 matWorld = glm::mat4(1.0f);
		glm::mat4 mvp = m_camera.GetViewProj() *matWorld;

		m_program_tess.SetUniform( "MVP", mvp );

		m_vbT.On();
		m_vbT.SetPatchVertices(2);
		m_vbT.Draw(GL_PATCHES, 0, 2*(N+1));

		m_vbT.Off();

		m_program_tess.Off();
	}
	else
	{
		m_program_curve.On();

		glm::mat4 matWorld = glm::mat4(1.0f);
		glm::mat4 mvp = m_camera.GetViewProj() *matWorld;

		m_program_curve.SetUniform( "MVP", mvp );

		m_vb.On();

		m_vb.Draw(GL_LINE_STRIP, 0, N+1);

		m_vb.Off();

		m_program_curve.Off();
	}
	m_program_basic.On();

	glm::mat4 matWorld = glm::mat4(1.0f);
	glm::mat4 mvp = m_camera.GetViewProj() *matWorld;

	m_program_basic.SetUniform( "MVP", mvp );

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
		std::cout << "[app.Init] Az alkalmazás inicializálása közben hibatörtént!" << std::endl;
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
		std::cout << "[SDL indítása]Hiba az SDL inicializálása közben: " << SDL_GetError() << std::endl;
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
		std::cout << "[Ablak létrehozása]Hiba az SDL inicializálása közben: " << SDL_GetError() << std::endl;
		exit(1);
	}

	context = SDL_GL_CreateContext(win);
	if (context == 0)
	{
		std::cout << "[OGL context létrehozása]Hiba az SDL inicializálása közben: " << SDL_GetError() << std::endl;
		exit(1);
	}

	SDL_GL_SetSwapInterval(1);

	GLenum error = glewInit();
	if (error != GLEW_OK)
	{
		std::cout << "[GLEW] Hiba az inicializálás során!" << std::endl;
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

		std::cout << "[OGL context létrehozása] Nem sikerült létrehozni az OpenGL context-et! Lehet, hogy az SDL_GL_SetAttribute(...) hívásoknál az egyik beállítás helytelen." << std::endl;

		exit(1);
	}

	std::stringstream window_title;
	window_title << "OpenGL " << glVersion[0] << "." << glVersion[1];
	SDL_SetWindowTitle(win, window_title.str().c_str());

}
