#pragma once

#include "../CurveLib/bezierCurve.h"
#include "../CurveLib/Curve.hpp"
#include "../CurveLib/SimpleCurves.hpp"

#include <vector>

class ExampleHandler
{
	static std::vector<Curve*> examples;
	
	static bool ready;

public:
	static Curve& get(int i);
};
