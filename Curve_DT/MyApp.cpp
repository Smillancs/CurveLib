#include "MyApp.h"
#include "GLUtils.hpp"

#include <glm/gtc/type_ptr.hpp>
//#include <GL/GLU.h>
//#include <math.h>
#include <sstream>

CMyApp::CMyApp()
{
	win = 0;
}

CMyApp::~CMyApp(void)
{
	delete currentCurve;
	delete generator;


}


bool CMyApp::Init()
{
	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
		
	generator = new CurveRenderer(*currentCurve);

	generator->genBufferCps(static_cast<BezierCurve*>(currentCurve)->GetGlmControlPoints());
	m_vbBez = generator->getBuffer();

	m_vbQuad.AddAttribute(0, 3);

	m_vbQuad.AddData(0, -100, 100, 0);
	m_vbQuad.AddData(0, 100, -100, 0);
	m_vbQuad.AddData(0, -100, -100, 0);
	m_vbQuad.AddData(0, 100, 100, 0);
		
	m_vbQuad.AddIndex(2, 3, 0);
	m_vbQuad.AddIndex(2, 1, 3);

	res.x = res.y = 20;

	m_vbQuad.InitBuffers();

	m_program_bez.AttachShader(GL_VERTEX_SHADER, "../Assets/shader_tess_adv.vert");
	m_program_bez.AttachShader(GL_TESS_CONTROL_SHADER, "../Assets/shader_tess_adv.tcs");
	m_program_bez.AttachShader(GL_TESS_EVALUATION_SHADER, "../Assets/shader_tess_adv.tes");
	m_program_bez.AttachShader(GL_FRAGMENT_SHADER, "../Assets/shader_tess_adv.frag");

	if ( !m_program_bez.LinkProgram())
	{
		return false;
	}

	m_camera.SetProj(45.0f, 640.0f/480.0f, 0.01f, 1000.0f);
	m_camera.SetView(glm::vec3(0,0,15), glm::vec3(0,0,0), glm::vec3(0,1,0));

	return true;
}

void CMyApp::Clean()
{
	delete generator;
	m_program_bez.Clean();
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

	auto cp = static_cast<BezierCurve*>(currentCurve)->GetGlmControlPoints();
	int size = cp.size();

	auto & currentProgram = m_program_bez;
	currentProgram.On();
	currentProgram.SetUniform("VP", vp);
	currentProgram.SetUniform("pointNum", size);

	glUniform3fv(glGetUniformLocation(currentProgram.GetProgramId(), "pointData"),
		3*cp.size(),
		&cp[0][0]);

	m_vbBez.On();
	m_vbBez.SetPatchVertices(vertsInPatch);
	m_vbBez.Draw(GL_PATCHES, 0, vertsInPatch);
	m_vbBez.Off();

	currentProgram.Off();
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
