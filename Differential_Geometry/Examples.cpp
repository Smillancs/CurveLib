#include "Examples.hpp"

bool ExampleHandler::ready = false;
std::vector<Curve*> ExampleHandler::examples;


Curve& ExampleHandler::get(int i)
{
	if(!ready)
	{
		examples.push_back(new Line(glm::dvec3(0,0,0), glm::dvec3(1,2,3)));
		examples.push_back(new Circle(glm::dvec3(1,2,3), 3.0));
		examples.push_back(new Spiral(glm::dvec3(0, 0, 0), 5));

		std::vector<ControlPoint> cps;
		cps.push_back({0,0,10});
		cps.push_back({0,5,0});
		cps.push_back({0,10,0});
		cps.push_back({-10,0,0});

		examples.push_back(new BezierCurve(cps));
		ready = true;
	}
	return *examples[i];
}
