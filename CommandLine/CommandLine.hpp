#pragma once

#include "../CurveLib/Utils.hpp"
#include "../CurveLib/Examples.hpp"
#include "../CurveLib/GeomInvariant.hpp"
#include "../CurveLib/SimpleCurves.hpp"
#include "../CurveLib/bezierCurve.h"

#include <iostream>
#include <vector>
#include <sstream>

class CommandLine
{
public:
	void run();

private:
	bool process(const std::string& cmd);
	
	bool curveProcess(std::stringstream& cmd);
	
	bool runTest();
	
	bool printHelp();
	
	size_t activeCurve = 0;
};

void CommandLine::run()
{
	bool ok = true;
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
		return runTest();
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
	else if(cmd1 == "T")
	{
		cmd >> cmd1;
		const char* s = cmd1.c_str();
		double t;
		sscanf(s, "%lf", &t);
		std::cout << GeomInv::T(ExampleHandler::get(activeCurve), t) << std::endl;
	}
	else if(cmd1 == "help")
	{
		std::cout << "Operations on curves:\n\tset - set active curve (from examples)\n\tinfo - print information about active curve\n\tf - evaluate curve in given parameter\n\td, dd, ddd - evaluate derivative of active curve\n\te, n, b - evaluate Frenet-frame of active curve\n\tK, T - evaluate curvature and torsion of active curve" << std::endl;
	}
	else
		throw Exception("This command does not exist for curves");
	return true;
}


bool CommandLine::runTest()
{
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

	std::cout << "All tests ran in order." << std::endl;
	return true;
}


bool CommandLine::printHelp()
{
	std::cout << "Available commands: \n\tcurve - perform operations on curves (`curve help` for info)\n\ttest - run unit tests\n\tquit - quit program\n";
	return true;
}
