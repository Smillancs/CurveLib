#include "Examples.hpp"

#include "SimpleCurves.hpp"
#include "bezierCurve.h"


bool ExampleHandler::ready = false;
std::vector<Curve*> ExampleHandler::examples;


Curve& ExampleHandler::get(int i)
{
	if(!ready)
	{
		generate();
	}
	return *examples[i];
}

size_t ExampleHandler::size() 
{
	if(!ready)
	{
		generate();
	}
	return examples.size(); 
}

void ExampleHandler::generate()
{
//	examples.push_back(new Line(glm::dvec3(0,0,0), glm::dvec3(1,2,3)));
//	examples.push_back(new Circle(glm::dvec3(1,2,3), 3.0));
//	examples.push_back(new Spiral(glm::dvec3(0, 0, 0), 5));

	std::vector<glm::vec3> cps;
	cps.push_back(glm::vec3(2,8,1));
	cps.push_back(glm::vec3(4,0,1));
	cps.push_back(glm::vec3(0,8,1));
	cps.push_back(glm::vec3(-5,0,1));
	cps.push_back(glm::vec3(1,-4,1));
	
	examples.push_back(new BezierCurve(cps));

	std::vector<glm::vec3> cps2;
	cps2.push_back(glm::vec3(1, 0, 1));
	cps2.push_back(glm::vec3(2, 3, 1));
	cps2.push_back(glm::vec3(-1, 1, 1));
	cps2.push_back(glm::vec3(0, -2, 1));
	cps2.push_back(glm::vec3(0, -2, 1));

	examples.push_back(new BezierCurve(cps2));

	std::vector<glm::vec3> cps3;
	cps3.push_back(glm::vec3(-1, -1, 1));
	cps3.push_back(glm::vec3(-1, 1, 1));
	cps3.push_back(glm::vec3(1, 1, 1));
	cps3.push_back(glm::vec3(1, -1, 1));
	cps3.push_back(glm::vec3(-1, -1, 1));

	examples.push_back(new BezierCurve(cps3));

	std::vector<glm::vec3> cps4;
	cps4.push_back(glm::vec3(-4, -4, 1));
	cps4.push_back(glm::vec3(-4, 4, 1));
	cps4.push_back(glm::vec3(4, 4, 1));
	cps4.push_back(glm::vec3(4, -4, 1));
	cps4.push_back(glm::vec3(8, -4, 1));

	examples.push_back(new BezierCurve(cps4));
	ready = true;
}
