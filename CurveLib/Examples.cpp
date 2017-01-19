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
	cps.push_back(glm::vec3(0,0,10));
	cps.push_back(glm::vec3(0,5,0));
	cps.push_back(glm::vec3(0,10,0));
	cps.push_back(glm::vec3(-10,0,0));

	examples.push_back(new BezierCurve(cps));
	ready = true;
}
