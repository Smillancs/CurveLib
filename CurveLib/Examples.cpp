#include "Examples.hpp"

#include "SimpleCurves.hpp"
#include "bezierCurve.h"
#include "RandomCurve.hpp"

bool ExampleHandler::ready = false;
std::vector<Curve::Ptr> ExampleHandler::examples;
Curve::Ptr ExampleHandler::random;

Curve& ExampleHandler::get(int i)
{
	if(!ready)
	{
		generate();
	}
  if(i == -1)
  {
    return *random;
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
	examples.push_back(Curve::Ptr(new Line(glm::dvec3(0,0,0), glm::dvec3(1,2,3))));
	examples.push_back(Curve::Ptr(new Circle(glm::dvec3(1,2,3), 3.0)));
	examples.push_back(Curve::Ptr(new Spiral(glm::dvec3(0, 0, 0), 5)));

	std::vector<glm::vec3> cps;
	cps.push_back(glm::vec3(0, 0, 4));
	cps.push_back(glm::vec3(0, 4, 0));
	cps.push_back(glm::vec3(-8, 2, 0));
	cps.push_back(glm::vec3(0, 4, -2));
	cps.push_back(glm::vec3(-2, -8, 4));
	cps.push_back(glm::vec3(0, 10, 0));

	examples.push_back(Curve::Ptr(new BezierCurve(cps)));

	std::vector<glm::vec3> cps2;
	cps2.push_back(glm::vec3(1,0,0));
	cps2.push_back(glm::vec3(2,3,1));
	cps2.push_back(glm::vec3(-1,1,2));
	cps2.push_back(glm::vec3(0,-2,0));
	cps2.push_back(glm::vec3(0,-2,0));

	examples.push_back(Curve::Ptr(new BezierCurve(cps2)));

  random = RandomCurve(3,3);
	ready = true;
}

void ExampleHandler::newRandom(int deg, int dim)
{
  random = RandomCurve(deg, dim);
}
