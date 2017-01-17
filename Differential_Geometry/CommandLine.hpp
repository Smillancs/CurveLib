#pragma once

#include "../CurveLib/Utils.hpp"
#include "Examples.hpp"
#include "../CurveLib/GeomInvariant.hpp"
#include "../CurveLib/SimpleCurves.hpp"

#include <iostream>
#include <vector>

class CommandLine
{
public:
	void run();

private:
	bool process(const std::string& cmd);
	
	bool curveProcess(const std::string& cmd);
	
	bool runTest();
	
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
	if(cmd.length()>6 && cmd.substr(0,6)==std::string("curve "))
		return curveProcess(cmd.substr(6));
	if(cmd=="test")
		return runTest();
	if(cmd=="quit")
		return false;
	else
		throw Exception("This command cannot be used");
}

bool CommandLine::curveProcess(const std::string& cmd)
{
	if(cmd.length()>3 && cmd.substr(0,3)==std::string("set"))
	{
		activeCurve = atoi(cmd.substr(4).c_str());
		std::cout << activeCurve << std::endl;
	}
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
	assert_double_equal(K(c,0.6), 3.0);
	assert_double_equal(T(c,0.4), 0.0);
	
	std::cout << "All tests ran in order." << std::endl;
	return true;
}
