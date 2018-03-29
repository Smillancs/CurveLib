#pragma once

#include "../CurveLib/Utils.hpp"
#include "../CurveLib/Examples.hpp"
#include "../CurveLib/GeomInvariant.hpp"
#include "../CurveLib/SimpleCurves.hpp"
#include "../CurveLib/bezierCurve.h"
#include "../CurveLib/RandomCurve.hpp"

#include "../GPUcompute/GeomOptimize.hpp"
#include "../GPUcompute/Reconstruction.hpp"

#include <SDL.h>
#include <SDL_opengl.h>
#include <GL/glew.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <chrono>
#include <random>

#include <glm/gtc/matrix_access.hpp>

#ifdef _WIN32
#define sscanf sscanf_s
#endif

class CommandLine
{
public:
	void run();

private:
	bool process(const std::string& cmd);

	bool curveProcess(std::stringstream& cmd);

	bool runTest(std::stringstream& cmd);

	void doStuff();

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
	if(cmd1 == "w"){
		doStuff();
		return true;
	}
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
    if(target=="") target = "curvature";

		/*std::vector<GeomOptimize::Input2D3> vec = {{glm::vec2(0,0),glm::vec2(1,0),1,-1}};
		GeomOptimize opt;
		std::shared_ptr<std::vector<float>> dump = std::shared_ptr<std::vector<float>>(new std::vector<float>(100));

		auto start = std::chrono::high_resolution_clock::now();
		std::vector<GeomOptimize::Result> res = opt.optimize2D3(target, vec, dump);
		auto end = std::chrono::high_resolution_clock::now();
		Curve::Ptr optCurve = opt.createResultCurve(vec[0], res[0]);

		std::cerr << "Single optimization by " << target << " done in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< " milliseconds, with results (" << res[0].t0 << "; " << res[0].t1 << ") generating a curve with norm " << res[0].norm << std::endl;

		assert_double_equal((double)(*dump)[0], 42.0); // dummy assert for buffer binding problems
		assert_(target != "curvatureD" || res[0].norm <= 1.0f); // Depends on initial configuration, but with the current it can be done.*/



		std::vector<Reconstruction<2>::Input> vec = {getPointData2(ExampleHandler::getP(3),0), getPointData2(ExampleHandler::getP(3),1)};
		Reconstruction<2> opt(false);
		std::shared_ptr<std::vector<float>> dump = std::shared_ptr<std::vector<float>>(new std::vector<float>(100));

		auto start = std::chrono::high_resolution_clock::now();
		std::vector<Reconstruction<2>::Result> res = opt.optimize(target, vec, dump);
		auto end = std::chrono::high_resolution_clock::now();
		Curve::Ptr optCurve = opt.createResultCurve(res[0]);

    /*std::cerr << "Debug: " << std::endl;
      for(int i=0;i<12;++i) std::cerr << (*dump)[i] << std::endl;*/

		std::cerr << "Single optimization by " << target << " done in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
			<< " milliseconds, with resulting curve: " << optCurve->about() << " , with norm " << res[0].second[0] << std::endl;

		assert_double_equal((double)(*dump)[0], 42.0); // dummy assert for buffer binding problems
		//assert_(target != "curvatureD" || res[0].norm <= 1.0f); // Depends on initial configuration, but with the current it can be done.
	}
	else if("help" == option)
	{
		std::cout << "Running tests:\n\t<no option> - run basic assertion tests on curve functions\n\topt - also test curve optimizations" << std::endl;
	}
	std::cout << "All tests ran in order." << std::endl;
	return true;
}

void CommandLine::doStuff()
{
	// generating grids (2d)
	/*std::ofstream file("curv4d10.txt");

	std::ofstream file("grid_bounded.txt");
	const int N = 20;
	std::vector<Reconstruction<1>::Input> vec;
	for(int i=0;i<=N;++i) for(int j=0;j<=N;++j)
	{
		// depends on interval of parameters
		float alpha = M_PI*i/N - M_PI_2;//2*M_PI*(i)/N;
		float beta = M_PI*j/N - M_PI_2;//2*M_PI*(j)/N;
		glm::vec3 p0 = glm::vec3(0,0,0);
		glm::vec3 e0 = glm::vec3(cos(alpha), sin(alpha), 0);
		glm::vec3 p1 = glm::vec3(1,0,0);
		glm::vec3 e1 = glm::vec3(cos(beta), sin(beta), 0);
		vec.push_back(ReconstructionData<1>(p0, e0));
		vec.push_back(ReconstructionData<1>(p1, e1));
	}
	Reconstruction<1,0> opt(false);
	std::function<double(Curve&,double)> optTarget = std::function<double(Curve&,double)>([](Curve& c, double t){return abs(glm::dot(c.dnf(t,1), c.dnf(t,2)));});

	std::string _optTarget = "const_velocity";

	// CPU/GPU
	std::vector<Reconstruction<1,0>::Result_cpu> res;
	//std::vector<Reconstruction<1,0>::Result> res;
	
	// CPU/GPU
	res = opt.optimize_cpu_alglib(optTarget, vec, nullptr);
	//res = opt.optimize(_optTarget, vec, nullptr);

	for(int i=0;i<=N;++i) for(int j=0;j<=N;++j)
	{
		float alpha = M_PI*i/N - M_PI_2;//2*M_PI*(i)/N;
		float beta = M_PI*j/N - M_PI_2;//2*M_PI*(j)/N;
		
		Curve::Ptr c = opt.createResultCurve(res[i*(N+1)+j]);

		file << alpha << " " << beta << " " << GeomInv::frenet_coordinate<1,1>(*c, 0) << " " << GeomInv::frenet_coordinate<1,1>(*c, 1) << " " << res[i*(N+1)+j].second << std::endl;
	}*/

	// generating grids (2d, only curvature, fixed angles)
	/*std::ofstream file("curv20.txt");
	for(int i=0;i<N;++i) for(int j=0;j<N;++j)
	{
		float alpha = M_PI / 4;
		float beta = -M_PI / 4;
		glm::vec3 p0 = glm::vec3(0,0,0);
		glm::vec3 e0 = glm::vec3(cos(alpha), sin(alpha), 0);
		glm::vec3 n0 = glm::vec3(cos(alpha-M_PI/2), sin(alpha-M_PI/2), 0);
		glm::vec3 p1 = glm::vec3(1,0,0);
		glm::vec3 e1 = glm::vec3(cos(beta), sin(beta), 0);
		glm::vec3 n1 = glm::vec3(cos(beta-M_PI/2), sin(beta-M_PI/2), 0);

		float k0 = i*2.0/N;
		float k1 = j*2.0/N;

		std::vector<Reconstruction<2>::Input> vec = {ReconstructionData<2>(p0, e0, n0, k0), ReconstructionData<2>(p1, e1, n1, k1)};

		Reconstruction<2,0> opt(false);
		std::function<double(Curve&,double)> optTarget = std::function<double(Curve&,double)>([](Curve& c, double t){return abs(glm::dot(c.dnf(t,1), c.dnf(t,2)));});

		std::vector<Reconstruction<2,0>::Result_cpu> res;
		
		res = opt.optimize_cpu_alglib(optTarget, vec, nullptr);
		Curve::Ptr c = opt.createResultCurve(res[0]);

		file << k0 << " " << k1 << " " << GeomInv::frenet_coordinate<1,1>(*c, 0) << " " << GeomInv::frenet_coordinate<2,1>(*c, 0) << " " << GeomInv::frenet_coordinate<1,1>(*c, 1) << " " << GeomInv::frenet_coordinate<2,1>(*c, 1) << " " << res[0].second << std::endl;

	}*/

	// generating grids (4d)
	/*std::ofstream file("curv4d20.txt");
	const int N = 10, M = 20;
	for(int i=0;i<=N;++i) for(int j=0;j<=N;++j) for(int k=-M/2;k<=M/2;++k) for(int l=-M/2;l<=M/2;++l)
	{
		float alpha = -M_PI/2+M_PI*(i)/N;
		float beta = -M_PI/2+M_PI*(j)/N;
		glm::vec3 p0 = glm::vec3(0,0,0);
		glm::vec3 e0 = glm::vec3(cos(alpha), sin(alpha), 0);
		glm::vec3 n0 = glm::vec3(cos(alpha-M_PI/2), sin(alpha-M_PI/2), 0);
		glm::vec3 p1 = glm::vec3(1,0,0);
		glm::vec3 e1 = glm::vec3(cos(beta), sin(beta), 0);
		glm::vec3 n1 = glm::vec3(cos(beta-M_PI/2), sin(beta-M_PI/2), 0);

		float k0 = k*4.0/M;
		float k1 = l*4.0/M;

		std::vector<Reconstruction<2>::Input> vec = {ReconstructionData<2>(p0, e0, n0, k0), ReconstructionData<2>(p1, e1, n1, k1)};

		Reconstruction<2,0> opt(false);
		std::function<double(Curve&,double)> optTarget = std::function<double(Curve&,double)>([](Curve& c, double t){return abs(glm::dot(c.dnf(t,1), c.dnf(t,2)));});

		std::vector<Reconstruction<2,0>::Result_cpu> res;
		
		res = opt.optimize_cpu_alglib(optTarget, vec, nullptr);
		Curve::Ptr c = opt.createResultCurve(res[0]);

		file << alpha << " " << beta << " " << k0 << " " << k1 << " " << GeomInv::frenet_coordinate<1,1>(*c, 0) << " " << GeomInv::frenet_coordinate<2,1>(*c, 0) << " " << GeomInv::frenet_coordinate<1,1>(*c, 1) << " " << GeomInv::frenet_coordinate<2,1>(*c, 1) << " " << res[0].second << std::endl;

	}*/

	// generating scattered data
	/*unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator (seed);
	std::ofstream scattered4d("scattered4d.txt");
	const int N = 100;
	for(int i=0;i<=N;++i)
	{
		std::uniform_real_distribution<double> angles(-M_PI_2, M_PI_2);
		std::uniform_real_distribution<double> curvs(-2, 2);
		float alpha = angles(generator);
		float beta = angles(generator);
		glm::vec3 p0 = glm::vec3(0,0,0);
		glm::vec3 e0 = glm::vec3(cos(alpha), sin(alpha), 0);
		glm::vec3 n0 = glm::vec3(cos(alpha-M_PI/2), sin(alpha-M_PI/2), 0);
		glm::vec3 p1 = glm::vec3(1,0,0);
		glm::vec3 e1 = glm::vec3(cos(beta), sin(beta), 0);
		glm::vec3 n1 = glm::vec3(cos(beta-M_PI/2), sin(beta-M_PI/2), 0);

		float k0 = curvs(generator);
		float k1 = curvs(generator);

		std::vector<Reconstruction<2>::Input> vec = {ReconstructionData<2>(p0, e0, n0, k0), ReconstructionData<2>(p1, e1, n1, k1)};

		Reconstruction<2,0> opt(false);
		std::function<double(Curve&,double)> optTarget = std::function<double(Curve&,double)>([](Curve& c, double t){return abs(glm::dot(c.dnf(t,1), c.dnf(t,2)));});

		std::vector<Reconstruction<2,0>::Result_cpu> res;
		
		res = opt.optimize_cpu_alglib(optTarget, vec, nullptr);
		Curve::Ptr c = opt.createResultCurve(res[0]);

		scattered4d << alpha << " " << beta << " " << k0 << " " << k1 << " " << GeomInv::frenet_coordinate<1,1>(*c, 0) << " " << GeomInv::frenet_coordinate<2,1>(*c, 0) << " " << GeomInv::frenet_coordinate<1,1>(*c, 1) << " " << GeomInv::frenet_coordinate<2,1>(*c, 1) << " " << res[0].second << std::endl;

	}*/


	// testing estimates by method
	/*std::vector<std::string> infiles = {"estimates.txt", "estimates_nn.txt", "estimates_cs.txt", "estimates_pch.txt", "estimates_kappa_pch.txt"};
	std::vector<std::string> outfiles = {"error.txt", "error_nn.txt", "error_cs.txt", "error_pch.txt", "error_kappa_pch.txt"};
	std::vector<std::string> exacts = {"exact.txt", "exact_nn.txt", "exact_cs.txt", "exact_pch.txt", "exact_kappa_pch.txt"};
	int i = 2;
	//std::cin >> i;
	std::ifstream est(infiles[i].c_str());
	std::ofstream err(outfiles[i].c_str());
	std::ofstream ex(exacts[i].c_str());

	for(int i=0;i<10000;++i)
	{
		float alpha, beta, x1b, x1j, wEst;
		est >> alpha >> beta >> x1b >> x1j >> wEst;
		glm::vec3 p0 = glm::vec3(0,0,0);
		glm::vec3 e0 = glm::vec3(cos(alpha), sin(alpha), 0);
		glm::vec3 p1 = glm::vec3(1,0,0);
		glm::vec3 e1 = glm::vec3(cos(beta), sin(beta), 0);
		std::vector<Reconstruction<1>::Input> vec = {ReconstructionData<1>(p0, e0), ReconstructionData<1>(p1, e1)};

		Reconstruction<1,0> opt(false);
		std::function<double(Curve&,double)> optTarget = std::function<double(Curve&,double)>([](Curve& c, double t){return abs(glm::dot(c.dnf(t,1), c.dnf(t,2)));});// GeomInv::dK;
		//std::function<double(Curve&,double)> optTarget = GeomInv::K;

		std::vector<Reconstruction<1,0>::Result_cpu> res;
		
		res = opt.optimize_cpu_alglib(optTarget, vec, nullptr);
		Curve::Ptr c = opt.createResultCurve(res[0]);

		//std::cout << alpha << " " << beta << " " << GeomInv::frenet_coordinate<1,1>(*c, 0) << " " << GeomInv::frenet_coordinate<1,1>(*c, 1) << " " << integrate(optTarget, c) << std::endl;
		//std::cout << c->about() << std::endl;

		float x1bOpt = GeomInv::frenet_coordinate<1,1>(*c, 0), x1jOpt = GeomInv::frenet_coordinate<1,1>(*c, 1);
		float wOpt = integrate(optTarget, c);
		//float x1bOpt, x1jOpt, wOpt;
		//ex >> alpha >> beta >> x1bOpt >> x1jOpt >> wOpt;

		//std::cerr << c->about() << std::endl;


		std::array<float,2> x1s = {x1b, x1j};
		std::array<glm::vec3,4> points = calculateControlPoints<1,0,4>(PointDerivatives<1,0,2>(vec, x1s, 0, false), PointDerivatives<1,0,2>(vec, x1s, 0, true), std::array<glm::vec3,0>());
		c = Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end())));

		float wReal = integrate(optTarget, c);

		//if(x1bOpt > 1e-6 && x1jOpt > 1e-6 && wReal > 1e-6 && wOpt > 1e-6)
		ex << alpha << ' ' << beta << ' ' << x1bOpt << ' ' << x1jOpt << ' ' << wOpt << std::endl;

		err << (x1b-x1bOpt)/std::max(abs(x1bOpt), abs(x1b)) << ' ' << (x1j-x1jOpt)/std::max(abs(x1jOpt), abs(x1j)) << ' ' << (wEst-wReal)/std::max(abs(wReal), abs(wEst)) << ' ' << (wReal-wOpt)/std::max(abs(wOpt), abs(wReal)) << std::endl;
		//std::cerr << alpha << " " << beta << " " << GeomInv::frenet_coordinate<1,1>(*c, 0) << " " << GeomInv::frenet_coordinate<1,1>(*c, 1) << " " << integrate(optTarget, c) << " " << wOpt << std::endl;
		//std::cout << c->about() << std::endl;
	}*/
	
	// testing estimates by gridsize
	/*std::vector<std::string> infiles = {"estimates5.txt", "estimates10.txt", "estimates15.txt", "estimates20.txt", "estimates30.txt"};
	std::vector<std::string> outfiles = {"error5.txt", "error10.txt", "error15.txt", "error20.txt", "error30.txt"};
	//int i = 2;
	//std::cin >> i;
	for(int k=0;k<5;++k){
		std::ifstream est(infiles[k].c_str());
		std::ofstream err(outfiles[k].c_str());

		for(int i=0;i<1000;++i)
		{
			float alpha, beta, x1b, x1j, wEst;
			est >> alpha >> beta >> x1b >> x1j >> wEst;
			glm::vec3 p0 = glm::vec3(0,0,0);
			glm::vec3 e0 = glm::vec3(cos(alpha), sin(alpha), 0);
			glm::vec3 p1 = glm::vec3(1,0,0);
			glm::vec3 e1 = glm::vec3(cos(beta), sin(beta), 0);
			std::vector<Reconstruction<1>::Input> vec = {ReconstructionData<1>(p0, e0), ReconstructionData<1>(p1, e1)};

			Reconstruction<1,0> opt(false);
			std::function<double(Curve&,double)> optTarget = std::function<double(Curve&,double)>([](Curve& c, double t){return abs(glm::dot(c.dnf(t,1), c.dnf(t,2)));});// GeomInv::dK;
			//std::function<double(Curve&,double)> optTarget = GeomInv::K;

			std::vector<Reconstruction<1,0>::Result_cpu> res;
			
			res = opt.optimize_cpu_alglib(optTarget, vec, nullptr);
			Curve::Ptr c = opt.createResultCurve(res[0]);

			//std::cout << alpha << " " << beta << " " << GeomInv::frenet_coordinate<1,1>(*c, 0) << " " << GeomInv::frenet_coordinate<1,1>(*c, 1) << " " << integrate(optTarget, c) << std::endl;
			//std::cout << c->about() << std::endl;

			float x1bOpt = GeomInv::frenet_coordinate<1,1>(*c, 0), x1jOpt = GeomInv::frenet_coordinate<1,1>(*c, 1);
			float wOpt = integrate(optTarget, c);
			//float x1bOpt, x1jOpt, wOpt;
			//ex >> alpha >> beta >> x1bOpt >> x1jOpt >> wOpt;

			//std::cerr << c->about() << std::endl;


			std::array<float,2> x1s = {x1b, x1j};
			std::array<glm::vec3,4> points = calculateControlPoints<1,0,4>(PointDerivatives<1,0,2>(vec, x1s, 0, false), PointDerivatives<1,0,2>(vec, x1s, 0, true), std::array<glm::vec3,0>());
			c = Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end())));

			float wReal = integrate(optTarget, c);

			//if(x1bOpt > 1e-6 && x1jOpt > 1e-6 && wReal > 1e-6 && wOpt > 1e-6)
			//ex << alpha << ' ' << beta << ' ' << x1bOpt << ' ' << x1jOpt << ' ' << wOpt << std::endl;

			err << (x1b-x1bOpt)/std::max(abs(x1bOpt), abs(x1b)) << ' ' << (x1j-x1jOpt)/std::max(abs(x1jOpt), abs(x1j)) << ' ' << (wEst-wReal)/std::max(abs(wReal), abs(wEst)) << ' ' << (wReal-wOpt)/std::max(abs(wOpt), abs(wReal)) << std::endl;
			//std::cerr << alpha << " " << beta << " " << GeomInv::frenet_coordinate<1,1>(*c, 0) << " " << GeomInv::frenet_coordinate<1,1>(*c, 1) << " " << integrate(optTarget, c) << " " << wOpt << std::endl;
			//std::cout << c->about() << std::endl;
		}
	}*/


	// testing estimates (4d)
	/*std::ifstream est("estimates4d.txt");
	std::ofstream err("error4d.txt");*/
	/*std::ifstream est("estimates_scattered.txt");
	std::ofstream err("error_scattered.txt");*/

	/*for(int i=0;i<10000;++i)
	{
		float alpha, beta, k0, k1, x1b, x2b, x1j, x2j, wEst;
		est >> alpha >> beta >> k0 >> k1 >> x1b >> x2b >> x1j >> x2j >> wEst;
		glm::vec3 p0 = glm::vec3(0,0,0);
		glm::vec3 e0 = glm::vec3(cos(alpha), sin(alpha), 0);
		glm::vec3 n0 = glm::vec3(cos(alpha-M_PI/2), sin(alpha-M_PI/2), 0);
		glm::vec3 p1 = glm::vec3(1,0,0);
		glm::vec3 e1 = glm::vec3(cos(beta), sin(beta), 0);
		glm::vec3 n1 = glm::vec3(cos(beta-M_PI/2), sin(beta-M_PI/2), 0);

		std::vector<Reconstruction<2>::Input> vec = {ReconstructionData<2>(p0, e0, n0, k0), ReconstructionData<2>(p1, e1, n1, k1)};

		Reconstruction<2,0> opt(false);
		std::function<double(Curve&,double)> optTarget = std::function<double(Curve&,double)>([](Curve& c, double t){return abs(glm::dot(c.dnf(t,1), c.dnf(t,2)));});// GeomInv::dK;

		std::vector<Reconstruction<2,0>::Result_cpu> res;
		
		res = opt.optimize_cpu_alglib(optTarget, vec, nullptr);
		Curve::Ptr c = opt.createResultCurve(res[0]);

		//std::cout << alpha << " " << beta << " " << GeomInv::frenet_coordinate<1,1>(*c, 0) << " " << GeomInv::frenet_coordinate<1,1>(*c, 1) << " " << integrate(optTarget, c) << std::endl;
		//std::cout << c->about() << std::endl;

		float x1bOpt = GeomInv::frenet_coordinate<1,1>(*c, 0), x1jOpt = GeomInv::frenet_coordinate<1,1>(*c, 1);
		float x2bOpt = GeomInv::frenet_coordinate<2,1>(*c, 0), x2jOpt = GeomInv::frenet_coordinate<2,1>(*c, 1);
		float wOpt = integrate(optTarget, c);

		//std::cerr << c->about() << std::endl;


		std::array<float,4> x1s = {x1b, x1j, x2b, x2j};
		std::array<glm::vec3,6> points = calculateControlPoints<2,0,6>(PointDerivatives<2,0,4>(vec, x1s, 0, false), PointDerivatives<2,0,4>(vec, x1s, 0, true), std::array<glm::vec3,0>());
		c = Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end())));

		float wReal = integrate(optTarget, c);

		if(x1bOpt > 1e-6 && x1jOpt > 1e-6 && abs(x2bOpt) > 1e-6 && abs(x2jOpt) > 1e-6 && abs(wReal) > 1e-6 && abs(wOpt) > 1e-6)
		err << ((x1b-x1bOpt)/std::max(abs(x1bOpt), abs(x1b))) << ' ' << (x2b-x2bOpt)/std::max(abs(x2bOpt), abs(x2b)) << ' ' << (x1j-x1jOpt)/std::max(abs(x1jOpt),abs(x1j)) << ' ' << (x2j-x2jOpt)/std::max(abs(x2jOpt),abs(x2j)) << ' ' << ((wEst-wReal)/std::max(abs(wReal), abs(wEst))) << ' ' << ((wReal-wOpt)/std::max(abs(wOpt), abs(wReal))) << std::endl;
		//std::cout << alpha << " " << beta << " " << GeomInv::frenet_coordinate<1,1>(*c, 0) << " " << GeomInv::frenet_coordinate<1,1>(*c, 1) << " " << integrate(optTarget, c) << std::endl;
		//std::cout << c->about() << std::endl;
	}*/
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
