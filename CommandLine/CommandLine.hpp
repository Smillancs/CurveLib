#pragma once

#include "../CurveLib/Utils.hpp"
#include "../CurveLib/Examples.hpp"
#include "../CurveLib/GeomInvariant.hpp"
#include "../CurveLib/SimpleCurves.hpp"
#include "../CurveLib/bezierCurve.h"

#include "../GPUcompute/GeomOptimize.hpp"

#include <SDL.h>
#include <SDL_opengl.h>
#include <GL/glew.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <chrono>

#include <glm/gtc/matrix_access.hpp>

class CommandLine
{
public:
	void run();

private:
	bool process(const std::string& cmd);

	bool curveProcess(std::stringstream& cmd);

	bool runTest(std::stringstream& cmd);

	bool printHelp();

	void initOGL();

	bool OGLinited = false;

	int activeCurve = 0;
};

void CommandLine::run()
{
	bool ok = true;
	initOGL();
	while(ok)
	{
		std::string cmd;
		try{
			std::cout << "> ";
			std::getline(std::cin, cmd);
			ok = process(cmd);
		}catch(std::exception& e){
			std::cout << "Error happened: " << e.what() << std::endl;
			ok = true;
		}
	}
}

bool CommandLine::process(const std::string& cmd)
{
	std::stringstream str ( cmd );
	std::string cmd1;
	str >> cmd1;
	if(cmd1 == "curve")
		return curveProcess(str);
	if(cmd1 == "test")
		return runTest(str);
	if(cmd1 == "help")
		return printHelp();
	if(cmd1 == "quit")
		return false;
	else
		throw Exception("This command cannot be used");
}

bool CommandLine::curveProcess(std::stringstream& cmd)
{
	std::string cmd1;
	cmd >> cmd1;
	if(cmd1 == "set")
	{
		cmd >> cmd1;
		int index = atoi(cmd1.c_str());
		if(index < ExampleHandler::size())
			activeCurve = index;
		std::cout << "Curve #" << activeCurve << ": " << ExampleHandler::get(activeCurve).about() << std::endl;
	}
  else if(cmd1 == "setrandom")
  {
    activeCurve = -1;
    std::cout << "Curve #-1 (random curve): " << ExampleHandler::get(activeCurve).about() << std::endl;
  }
  else if(cmd1 == "newrandom")
  {
    activeCurve = -1;
  	cmd >> cmd1;
    int deg = atoi(cmd1.c_str());
  	cmd >> cmd1;
    int dim = atoi(cmd1.c_str());
    ExampleHandler::newRandom(deg, dim);
    std::cout << "Curve #-1 (random curve): " << ExampleHandler::get(activeCurve).about() << std::endl;
  }
	else if(cmd1 == "info")
	{
		std::cout << "Curve #" << activeCurve << ": " << ExampleHandler::get(activeCurve).about() << std::endl;
	}
	else if(cmd1 == "f")
	{
		cmd >> cmd1;
		const char* s = cmd1.c_str();
		double t;
		sscanf(s, "%lf", &t);
		std::cout << glm::to_string(ExampleHandler::get(activeCurve).f(t)) << std::endl;
	}
	else if(cmd1 == "d")
	{
		cmd >> cmd1;
		const char* s = cmd1.c_str();
		double t;
		sscanf(s, "%lf", &t);
		std::cout << glm::to_string(ExampleHandler::get(activeCurve).dnf(t,1)) << std::endl;
	}
	else if(cmd1 == "dd")
	{
		cmd >> cmd1;
		const char* s = cmd1.c_str();
		double t;
		sscanf(s, "%lf", &t);
		std::cout << glm::to_string(ExampleHandler::get(activeCurve).dnf(t,2)) << std::endl;
	}
	else if(cmd1 == "ddd")
	{
		cmd >> cmd1;
		const char* s = cmd1.c_str();
		double t;
		sscanf(s, "%lf", &t);
		std::cout << glm::to_string(ExampleHandler::get(activeCurve).dnf(t,3)) << std::endl;
	}
	else if(cmd1 == "e")
	{
		cmd >> cmd1;
		const char* s = cmd1.c_str();
		double t;
		sscanf(s, "%lf", &t);
		std::cout << glm::to_string(GeomInv::e(ExampleHandler::get(activeCurve), t)) << std::endl;
	}
	else if(cmd1 == "n")
	{
		cmd >> cmd1;
		const char* s = cmd1.c_str();
		double t;
		sscanf(s, "%lf", &t);
		std::cout << glm::to_string(GeomInv::n(ExampleHandler::get(activeCurve), t)) << std::endl;
	}
	else if(cmd1 == "b")
	{
		cmd >> cmd1;
		const char* s = cmd1.c_str();
		double t;
		sscanf(s, "%lf", &t);
		std::cout << glm::to_string(GeomInv::b(ExampleHandler::get(activeCurve), t)) << std::endl;
	}
	else if(cmd1 == "K")
	{
		cmd >> cmd1;
		const char* s = cmd1.c_str();
		double t;
		sscanf(s, "%lf", &t);
		std::cout << GeomInv::K(ExampleHandler::get(activeCurve), t) << std::endl;
	}
	else if(cmd1 == "dK")
	{
		cmd >> cmd1;
		const char* s = cmd1.c_str();
		double t;
		sscanf(s, "%lf", &t);
		std::cout << GeomInv::dK(ExampleHandler::get(activeCurve), t) << std::endl;
	}
	else if(cmd1 == "ddK")
	{
		cmd >> cmd1;
		const char* s = cmd1.c_str();
		double t;
		sscanf(s, "%lf", &t);
		std::cout << GeomInv::ddK(ExampleHandler::get(activeCurve), t) << std::endl;
	}
	else if(cmd1 == "T")
	{
		cmd >> cmd1;
		const char* s = cmd1.c_str();
		double t;
		sscanf(s, "%lf", &t);
		std::cout << GeomInv::T(ExampleHandler::get(activeCurve), t) << std::endl;
	}
	else if(cmd1 == "dT")
	{
		cmd >> cmd1;
		const char* s = cmd1.c_str();
		double t;
		sscanf(s, "%lf", &t);
		std::cout << GeomInv::dT(ExampleHandler::get(activeCurve), t) << std::endl;
	}
	else if(cmd1 == "help")
	{
		std::cout << "Operations on curves:\n\tset - set active curve (from examples)\n\tsetrandom - set last generated random curve to active\n\tnewrandom <deg> <dim> - generate new random curve with given degree and dimension (2 or 3)\n\tinfo - print information about active curve\n\tf - evaluate curve in given parameter\n\td, dd, ddd - evaluate derivative of active curve\n\te, n, b - evaluate Frenet-frame of active curve\n\tK, dK, ddK - evaluate curvature and its derivatives\n\tT, dT - evaluate torsion and its derivatives" << std::endl;
	}
	else
		throw Exception("This command does not exist for curves");
	return true;
}


bool CommandLine::runTest(std::stringstream& cmd)
{
	std::string option;
	cmd >> option;

	using namespace GeomInv;
	Line l(glm::dvec3(0,0,0), glm::dvec3(1,1,0));
	assert_dvec3_equal(e(l,0), glm::dvec3(sqrt(0.5),sqrt(0.5),0));

	Circle c(glm::dvec3(1,2,3), 3.0);
	assert_dvec3_equal(e(c,0.5), glm::dvec3(0,-1,0));
	assert_dvec3_equal(n(c,0.5), glm::dvec3(1,0,0));
	assert_dvec3_equal(b(c,0.324), glm::dvec3(0,0,1));
	assert_double_equal(K(c,0.6), 0.333333);
	assert_double_equal(T(c,0.4), 0.0);

	//TODO: write tests for other bezier functions too
	std::vector<glm::vec3> cps;
	cps.push_back(glm::vec3(0, 2, 0));
	cps.push_back(glm::vec3(4, 0, 0));
	cps.push_back(glm::vec3(0, 8, 16));
	cps.push_back(glm::vec3(-4, 0, 2));
	BezierCurve bez(cps);

	auto first = bez.GetControlPoints().at(0);
	auto last = bez.GetControlPoints().at(bez.GetControlPoints().size()-1);

	assert_dvec3_equal(bez.f(0),glm::dvec3(first.x(),first.y(),first.z()));
	assert_dvec3_equal(bez.f(1),glm::dvec3(last.x(),last.y(),last.z()));
	assert_dvec3_equal(e(bez, 0.5), glm::dvec3(-0.388514, 0.291386, 0.874157));
	assert_dvec3_equal(n(bez, 0.5), glm::dvec3(-0.917284, -0.212334, -0.336903));
	assert_double_equal(T(bez, 0.5), 0.0860946);

	Curve& be = ExampleHandler::get(4);

	assert_dvec3_equal(e(be,0), glm::dvec3(0.301511, 0.904534, 0.301511));
	assert_dvec3_equal(n(be,0), glm::dvec3(-0.794552, 0.063564, 0.603860));
	assert_dvec3_equal(b(be,0), glm::dvec3(0.527046, -0.421637, 0.737865));
	assert_double_equal(K(be,0), 0.195026);
	assert_double_equal(T(be,0), 0.0166667);

	if("opt" == option)
	{
    std::string target;
    cmd >> target;

		std::vector<GeomOptimize::Input2D3> vec = {{glm::vec2(0,0),glm::vec2(1,0),1,-1}};
		GeomOptimize opt;
		std::shared_ptr<std::vector<float>> dump = std::shared_ptr<std::vector<float>>(new std::vector<float>(100));

		auto start = std::chrono::high_resolution_clock::now();
		std::vector<GeomOptimize::Result> res = opt.optimize2D3(target, vec, dump);
		auto end = std::chrono::high_resolution_clock::now();
		BezierCurve optCurve = opt.createResultCurve(vec[0], res[0]);

		std::cerr << "Single optimization by " << target << " done in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< " milliseconds, with results (" << res[0].t0 << "; " << res[0].t1 << ") generating a curve with norm " << res[0].norm << std::endl;

		assert_double_equal((double)(*dump)[0], 42.0); // dummy assert for buffer binding problems
		assert_(target != "curvatureD" || res[0].norm <= 1.0f); // Depends on initial configuration, but with the current it can be done.
	}
	else if("help" == option)
	{
		std::cout << "Running tests:\n\t<no option> - run basic assertion tests on curve functions\n\topt - also test curve optimizations" << std::endl;
	}
	std::cout << "All tests ran in order." << std::endl;
	return true;
}


bool CommandLine::printHelp()
{
	std::cout << "Available commands: \n\tcurve - perform operations on curves (`curve help` for info)\n\ttest - run unit tests\n\tquit - quit program\n";
	return true;
}

void CommandLine::initOGL()
{
	if(OGLinited) return;

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

	SDL_Window* win = SDL_CreateWindow("Hello SDL&OpenGL!",	100, 100, 640, 480,	SDL_WINDOW_OPENGL | SDL_WINDOW_HIDDEN | SDL_WINDOW_RESIZABLE);

	if (win == 0)
	{
		std::cout << "[Ablak létrehozása]Hiba az SDL inicializálása közben: " << SDL_GetError() << std::endl;
		exit(1);
	}

	SDL_GLContext context = SDL_GL_CreateContext(win);
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

	OGLinited = true;
}
