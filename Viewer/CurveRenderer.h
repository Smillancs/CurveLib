#pragma once

#include "../CurveLib/Curve.hpp"

#include "gVertexBuffer.h"

class CurveRenderer
{
public:
	CurveRenderer(Curve& curve) :c(curve){}

	void genBuffer(unsigned N, float a, float b);
	void changeCurve(Curve& curve){ c = curve; dirty = true; }

	gVertexBuffer getBuffer();

	Curve& c;

private:
	gVertexBuffer vb;
	unsigned n = 100;
	float a,b;
	bool dirty = true;
};
