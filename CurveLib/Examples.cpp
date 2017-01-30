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
	examples.push_back(new Line(glm::dvec3(0,0,0), glm::dvec3(1,2,3)));
	examples.push_back(new Circle(glm::dvec3(1,2,3), 3.0));
	examples.push_back(new Spiral(glm::dvec3(0, 0, 0), 5));

	std::vector<glm::vec3> cps;
	cps.push_back(glm::vec3(0, 0, 4));
	cps.push_back(glm::vec3(0, 4, 0));
	cps.push_back(glm::vec3(-8, 2, 0));
	cps.push_back(glm::vec3(0, 4, -2));
	cps.push_back(glm::vec3(-2, -8, 4));
	cps.push_back(glm::vec3(0, 10, 0));
	
	examples.push_back(new BezierCurve(cps));
	
	std::vector<glm::vec3> cps2;
	cps2.push_back(glm::vec3(1,0,0));
	cps2.push_back(glm::vec3(2,3,1));
	cps2.push_back(glm::vec3(-1,1,2));
	cps2.push_back(glm::vec3(0,-2,0));
	cps2.push_back(glm::vec3(0,-2,0));

	examples.push_back(new BezierCurve(cps2));
	ready = true;
}
