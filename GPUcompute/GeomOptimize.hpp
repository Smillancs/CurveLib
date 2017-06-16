#pragma once

#include "../CurveLib/Curve.hpp"
#include "../CurveLib/bezierCurve.h"

#include "gShaderProgram.h"

#include <vector>
#include <memory>

class GeomOptimize
{
public:
	struct Result
	{
		float t0;
		float t1;
		float norm;
	};

	struct Input2D3
	{
		glm::vec2 p0;
		glm::vec2 p1;
		float alpha;
		float beta;
	};

	std::vector<Result> optimize2D3(const std::string& targetFunction, const std::vector<Input2D3>& input, const std::shared_ptr<std::vector<float>>& debugInfo = 0);

	GeomOptimize();

	BezierCurve createResultCurve(const Input2D3& input, const Result& res);

protected:
	gShaderProgram program;
};
