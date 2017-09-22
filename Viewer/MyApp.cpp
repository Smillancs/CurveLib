#include "MyApp.h"
#include "GLUtils.hpp"

#include "../GPUcompute/GeomOptimize.hpp"
#include "../GPUcompute/Reconstruction.hpp"
#include "../CurveLib/RandomCurve.hpp"
#include "../CurveLib/Utils.hpp"

#include "../Imgui/imgui.h"
#include "imgui_impl_sdl_gl3.h"

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
  glLineWidth(2.0);

	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);

  axes = coordAxes();

	std::vector<GeomOptimize::Input2D3> vec = { {glm::vec2(0,0),glm::vec2(1,0),1,-1} };
	GeomOptimize opt;
	std::vector<GeomOptimize::Result> res = opt.optimize2D3("curvatureD", vec);
	optCurve = opt.createResultCurve(vec[0], res[0]);

  randomCurve = RandomCurve(5,3);

	/*Curve::Ptr c = optCurve;
	Curve::Ptr c = ExampleHandler::getP(3);
  Curve::Ptr c = randomCurve;*/

  InitCurveRenderer(randomCurve, true);
  activeCurve = randomCurve;
  /*InitCurveRenderer(ExampleHandler::getP(3), true);
  activeCurve = ExampleHandler::getP(3);*/

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

	m_camera.SetProj(45.0f, 4.0f/3.0f, 0.01f, 1000.0f);

	return true;
}

void CMyApp::InitCurveRenderer(Curve::Ptr c, bool gputesselation)
{
  m_vb.resize(1);
  m_vb[0].Clean();
  CurveRenderer ren(c);
  if(gputesselation)
  {
    vertsInPatch.resize(1);
    vertsInPatch[0] = ren.genBufferCps();
    m_vb[0] = ren.getBuffer();
    m_vb[0].InitBuffers();
    activeCurve = c;
  }
  else
  {
    ren.genBufferNormal(N, 0, 1);
    m_vb[0] = ren.getBuffer();
    m_vb[0].InitBuffers();
    activeCurve = c;
  }
}

void CMyApp::InitCurveRenderer(std::vector<Curve::Ptr> c, bool gputesselation)
{
  size_t n = c.size();
  m_vb.resize(n);
  vertsInPatch.resize(n);
  for(size_t i=0;i<n;++i)
  {
    m_vb[i].Clean();
    CurveRenderer ren(c[i]);
    if(gputesselation)
    {
      vertsInPatch[i] = ren.genBufferCps();
      m_vb[i] = ren.getBuffer();
      m_vb[i].InitBuffers();
    }
    else
    {
      ren.genBufferNormal(N, 0, 1);
      m_vb[i] = ren.getBuffer();
      m_vb[i].InitBuffers();
    }
  }
  activeCurves = c;
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

  if(change)
  {
    if(curveSelect == 0) activeCurve = optCurve;
    else activeCurve = ExampleHandler::getP(curveSelect);
    change = false;
    InitCurveRenderer(activeCurve, tesselated);
    activeCurves = {activeCurve};
  }

  if(makeOpt && optRank == 1)
  {
    std::vector<Reconstruction<1>::Input> vec = {getPointData1(activeCurve,0)};
    for(float i = 1; i <= segments; ++i)
    {
      vec.push_back(getPointData1(activeCurve,i/(float)segments));
      vec.push_back(getPointData1(activeCurve,i/(float)segments));
    }
    vec.pop_back();
		Reconstruction<1> opt;
		std::shared_ptr<std::vector<float>> dump = std::shared_ptr<std::vector<float>>(new std::vector<float>(100));

		std::vector<Reconstruction<1>::Result> res = opt.optimize("curvatureD", vec, dump);


        std::cerr << "Debug: " << std::endl;
          for(int i=0;i<12;++i) std::cerr << (*dump)[i] << std::endl;

    activeCurves.clear();
    for(size_t i=0;i<res.size();++i)
		{
      Curve::Ptr optCurve = opt.createResultCurve(res[i]);
      std::cerr << optCurve->about() << std::endl;
      activeCurves.push_back(optCurve);

      exportCurveData("curve.txt", optCurve, GeomInv::v, 1000, i);
    }


    makeOpt = false;

    InitCurveRenderer(activeCurves, tesselated);
  }
  if(makeOpt && optRank == 2)
  {
    std::vector<Reconstruction<2>::Input> vec = {getPointData2(activeCurve,0)};
    for(float i = 1; i <= segments; ++i)
    {
      vec.push_back(getPointData2(activeCurve,i/(float)segments));
      vec.push_back(getPointData2(activeCurve,i/(float)segments));
    }
    vec.pop_back();
		Reconstruction<2> opt;
		std::shared_ptr<std::vector<float>> dump = std::shared_ptr<std::vector<float>>(new std::vector<float>(100));

		std::vector<Reconstruction<2>::Result> res = opt.optimize("curvatureD", vec, dump);


        std::cerr << "Debug: " << std::endl;
          for(int i=0;i<12;++i) std::cerr << (*dump)[i] << std::endl;

    activeCurves.clear();
    for(size_t i=0;i<res.size();++i)
		{
      Curve::Ptr optCurve = opt.createResultCurve(res[i]);
      std::cerr << optCurve->about() << std::endl;
      activeCurves.push_back(optCurve);

      exportCurveData("curve.txt", optCurve, GeomInv::v, 1000, i);
    }


    makeOpt = false;

    InitCurveRenderer(activeCurves, tesselated);
  }
}


void CMyApp::Render()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glm::mat4 vp;
  if(ortho)
  {
    vp = glm::ortho(-zoom,zoom,-zoom,zoom,1.0f,20.0f)*glm::lookAt(glm::vec3(0,0,10),glm::vec3(0),glm::vec3(0,1,0));
  }
  else
    vp = m_camera.GetViewProj();

	if (tesselated)
	{
		auto & currentProgram = m_program_bez;
		currentProgram.On();
		currentProgram.SetUniform("VP", vp);
    for(size_t i=0;i<m_vb.size();++i)
    {
  		currentProgram.SetUniform("point_num", vertsInPatch[i]);
      std::shared_ptr<BezierCurve> bez = std::dynamic_pointer_cast<BezierCurve>(activeCurves.size() == 0 ? activeCurve : activeCurves[i]);
      auto cp = bez->GetGlmControlPoints();
      glUniform3fv(glGetUniformLocation(currentProgram.ID(), "pointData"),
      		cp.size(),
      		&cp[0][0]);
  		m_vb[i].On();
  		m_vb[i].SetPatchVertices(vertsInPatch[i]);
  		m_vb[i].Draw(GL_PATCHES, 0, vertsInPatch[i]);
  		m_vb[i].Off();
    }
		currentProgram.Off();
	}
	else
	{
		m_program_curve.On();

		m_program_curve.SetUniform( "VP", vp );

    for(size_t i=0;i<m_vb.size();++i)
    {
  		m_vb[i].On();
  		m_vb[i].Draw(GL_LINE_STRIP, 0, N+1);
  		m_vb[i].Off();
    }

		m_program_curve.Off();
	}
	m_program_basic.On();

	m_program_basic.SetUniform( "VP", vp );

	axes.On();
	axes.Draw(GL_LINES, 0, 6);
	axes.Off();

	m_program_basic.Off();


  ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiSetCond_FirstUseEver);
  ImGui::SetNextWindowSize(ImVec2(300,300), ImGuiSetCond_FirstUseEver);
	ImGui::Begin("Curve operations");
  if (ImGui::CollapsingHeader("Select curve"))
  {
      ImGui::TextWrapped("Select curve from prepared examples:");

      ImGui::RadioButton("Randomized curve: ", &curveSelect, -1);
      ImGui::RadioButton("Curve #1: ", &curveSelect, 3);
      ImGui::RadioButton("Curve #2: ", &curveSelect, 4);
      ImGui::RadioButton("Curve from simple optimization: ", &curveSelect, 0);

      if(ImGui::Button("Go!", ImVec2(0,0)))
      {
        change = true;
      }
  }
  if (ImGui::CollapsingHeader("Select drawing method"))
  {
      bool b = tesselated;
      ImGui::Checkbox("GPU tesselation", &tesselated);
      change |= (b != tesselated);
  }
  if (ImGui::CollapsingHeader("Camera settings"))
  {
      static int item = 0;
      ImGui::Combo("Projection", &item, "Free perspective\0Ortographic to XY\0\0");
      ortho = item == 1;
  }
  if (ImGui::CollapsingHeader("Generate new random curve"))
  {
      static char str0[5] = "5";
      ImGui::InputText("Degree", str0, 5);
      static int item = 1;
      ImGui::Combo("Dimension", &item, "2\0""3\0\0");
      if(ImGui::Button("Generate", ImVec2(0,0)))
      {
        ExampleHandler::newRandom(atoi(str0), item+2);
        change = true;
      }
  }
  if (ImGui::CollapsingHeader("Optimize current curve"))
  {
      static char str0[5] = "1";
      ImGui::InputText("Curve segments: ", str0, 5);
      segments = atoi(str0);
      if(ImGui::Button("Optimize (deg3)", ImVec2(0,0)))
      {
        makeOpt = true;
        optRank = 1;
      }
      if(ImGui::Button("Optimize (deg5)", ImVec2(0,0)))
      {
        makeOpt = true;
        optRank = 2;
      }
  }
	ImGui::End();

	/*ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiSetCond_FirstUseEver);
	ImGui::ShowTestWindow();*/

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
  zoom -= wheel.y * 0.1f;
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
	unsigned t = SDL_GetTicks();
	while (!quit)
	{
		while (SDL_PollEvent(&ev))
		{
			ImGui_ImplSdlGL3_ProcessEvent(&ev);
			bool is_mouse_captured = ImGui::GetIO().WantCaptureMouse;
			bool is_keyboard_captured = ImGui::GetIO().WantCaptureKeyboard;
			switch (ev.type)
				{
				case SDL_QUIT:
					quit = true;
					break;
				case SDL_KEYDOWN:
					if (ev.key.keysym.sym == SDLK_ESCAPE)
						quit = true;
					if (!is_keyboard_captured)
						KeyboardDown(ev.key);
					break;
				case SDL_KEYUP:
					if (!is_keyboard_captured)
						KeyboardUp(ev.key);
					break;
				case SDL_MOUSEBUTTONDOWN:
					if (!is_mouse_captured)
						MouseDown(ev.button);
					break;
				case SDL_MOUSEBUTTONUP:
					if (!is_mouse_captured)
						MouseUp(ev.button);
					break;
				case SDL_MOUSEWHEEL:
					if (!is_mouse_captured)
						MouseWheel(ev.wheel);
					break;
				case SDL_MOUSEMOTION:
					if (!is_mouse_captured)
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

    ImGui_ImplSdlGL3_NewFrame(win);

    Update();
		Render();
		ImGui::Render();

		SDL_GL_SwapWindow(win);
    unsigned spf = SDL_GetTicks() - t;
		if (spf < 16 ) SDL_Delay(16-spf);
  	t = SDL_GetTicks();
	}

	Clean();

	ImGui_ImplSdlGL3_Shutdown();

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

	win = SDL_CreateWindow("Hello SDL&OpenGL!",	100, 100, 1024, 768,	SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);

	if (win == 0)
	{
		std::cout << "[Ablak l�trehoz�sa]Hiba az SDL inicializ�l�sa k�zben: " << SDL_GetError() << std::endl;
		exit(1);
	}

	ImGui_ImplSdlGL3_Init(win);

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
