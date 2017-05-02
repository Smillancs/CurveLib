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

	m_program_lines.AttachShader(GL_VERTEX_SHADER, "../Assets/shader_basic.vert");
	m_program_lines.AttachShader(GL_FRAGMENT_SHADER, "../Assets/shader_basic.frag");

	if (!m_program_lines.LinkProgram())
	{
		return false;
	}

	m_camera.SetProj(45.0f, 640.0f/480.0f, 0.01f, 1000.0f);
	m_camera.SetView(glm::vec3(0,0,0), 10, glm::vec3(0,1,0));

	//int width = 640, height = 480;
	/*
	glGenFramebuffers(1, &fb);
	glGenTextures(1, &colorBuffer);

	glBindFramebuffer(GL_FRAMEBUFFER, fb);
	glBindTexture(GL_TEXTURE_2D, colorBuffer);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, NULL);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, colorBuffer, 0);
	if (glGetError() != GL_NO_ERROR)
	{
		std::cout << "Error creating color attachment" << std::endl;
		char ch; std::cin >> ch;
		exit(1);
	}


	glGenRenderbuffers(1, &pointBuffer);
	glBindRenderbuffer(GL_RENDERBUFFER, pointBuffer);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width, height);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, pointBuffer);
	if (glGetError() != GL_NO_ERROR)
	{
		std::cout << "Error creating depth attachment" << std::endl;
		char ch; std::cin >> ch;
		exit(1);
	}

	GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if (status != GL_FRAMEBUFFER_COMPLETE)
	{
		std::cout << "Incomplete framebuffer (";
		switch (status){
		case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:
			std::cout << "GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT";
			break;
		case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT:
			std::cout << "GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT";
			break;
		case GL_FRAMEBUFFER_UNSUPPORTED:
			std::cout << "GL_FRAMEBUFFER_UNSUPPORTED";
			break;
		}
		std::cout << ")" << std::endl;
		char ch;
		std::cin >> ch;
		exit(1);
	}

	// -- Unbind framebuffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	*/
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

glm::vec3 CMyApp::QuerryAndShow(const std::vector<glm::vec3> & cp, const int size, const glm::vec3 &point, const glm::mat4 &trf)
{
//	glBindFramebuffer(GL_FRAMEBUFFER, fb);
//	glClearColor(0.25f, 0.5f, 0.75f, 1.0f);
//	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	auto & currentProgram = m_program_bez;
	currentProgram.On();
	currentProgram.SetUniform("VP", trf);
	currentProgram.SetUniform("pointNum", size);
	currentProgram.SetUniform("referencePoint", point);
	glUniform3fv(glGetUniformLocation(currentProgram.GetProgramId(), "pointData"),
		3*cp.size(),
		&cp[0][0]);

	m_vbBez.On();
	m_vbBez.SetPatchVertices(vertsInPatch);
	m_vbBez.Draw(GL_PATCHES, 0, vertsInPatch);
	m_vbBez.Off();
	currentProgram.Off();
	
	//glm::vec4 result;
	//glReadPixels(0, 0, 1, 1, GL_RGBA, GL_FLOAT, &result);
/*	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	glViewport(0, 0, 640, 480);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glBindFramebuffer(GL_READ_FRAMEBUFFER, fb); 
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0); 
	glBlitFramebuffer(0, 0, 640, 480, 0, 0, 640, 480, GL_COLOR_BUFFER_BIT, GL_LINEAR); 
	glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);*/
	//return glm::vec3(result.x, result.y, result.z);
	return glm::vec3(0,0,0);
}

void CMyApp::DrawHull(const gVertexBuffer & from, int size, const glm::mat4 &trf)
{
	auto & currentProgram = m_program_lines;
	currentProgram.On();
	
	m_vbBez.On();
	m_vbBez.Draw(GL_LINES, 0, size*2+1);
	m_vbBez.Off();
	currentProgram.Off();
}

void CMyApp::DrawGuess(const glm::vec3 &point, const glm::vec3 &guess)
{
	auto & currentProgram = m_program_lines;
	currentProgram.On();
	gVertexBuffer temp;
	temp.On();
	temp.AddAttribute(0, 3);
	temp.AddData(0, point);
	temp.Draw(GL_POINTS, 0, 2);
	temp.Off();
	currentProgram.Off();
}

void CMyApp::Render()
{

	glClearColor(0.5, 0.5, 0.5, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//draw bezier
	glm::mat4 vp = m_camera.GetViewProj();
	auto cp = static_cast<BezierCurve*>(currentCurve)->GetGlmControlPoints();
	int size = cp.size();
	glm::vec3 point(0,-5,5);
	if (!cp.empty())
	{
		m_program_lines.On();
		m_program_lines.SetUniform("VP", vp);
		m_program_lines.SetUniform("pointNum", size);
		m_program_lines.SetUniform("rectRes", res);
		m_program_lines.SetUniform("zoom", m_camera.GetDistance());

		glUniform3fv(glGetUniformLocation(m_program_lines.GetProgramId(), "pointData"),
			3*cp.size(),
			&cp[0][0]);
		m_vbQuad.On();
		m_vbQuad.DrawIndexed(GL_TRIANGLES, 0, 6);
		m_vbQuad.Off();
		m_program_lines.Off();


		auto guess = QuerryAndShow(cp, size, point, vp);
		DrawHull(m_vbBez,size,vp);
		DrawGuess(point,guess);
	}
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
