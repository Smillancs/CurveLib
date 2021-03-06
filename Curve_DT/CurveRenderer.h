#pragma once

#include "../CurveLib/Curve.hpp"

#include "gVertexBuffer.h"

class CurveRenderer
{
public:
	CurveRenderer(Curve& curve) :c(curve){}

	void genBufferNormal(const unsigned N, const float a, const float b);
	void genBufferTesselation(const unsigned N, const float a, const float b);
	void genBufferCps(const std::vector<glm::vec3>&);
	void changeCurve(Curve& curve){ c = curve; dirty = true; }

	gVertexBuffer getBuffer();

	Curve& c;

private:
	gVertexBuffer vb;
	unsigned n = 100;
	float a,b;
	bool dirty = true;
};
