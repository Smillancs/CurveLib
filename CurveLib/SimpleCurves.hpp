#pragma once

#include "Curve.hpp"
#include "Utils.hpp"

// These curves are mostly for testing simple functions
// TODO: move all computation to dnf(t,n), and remove the dddf-like functions

class Line : public Curve
{
public:
	Line(glm::dvec3 a, glm::dvec3 b) :x(a), y(b){}

	glm::dvec3 dnf(double t, unsigned n){ switch (n){ case 1: return df(t); case 2: return ddf(t); case 3: return dddf(t); default: throw UnsupportedDerivativeOrder(n); } }

	glm::dvec3 f(double t){ return t*y + (1 - t)*x; }
	glm::dvec3 df(double t){ return glm::normalize(y - x); }
	glm::dvec3 ddf(double t){ return glm::dvec3(0.0); }
	glm::dvec3 dddf(double t){ return glm::dvec3(0.0); }

private:
	glm::dvec3 x;
	glm::dvec3 y;
};

class Circle : public Curve
{
public:
	Circle(glm::dvec3 p, double x) :o(p), r(x){}

	glm::dvec3 dnf(double t, unsigned n){ switch (n){ case 1: return df(t); case 2: return ddf(t); case 3: return dddf(t); default: throw UnsupportedDerivativeOrder(n); } }

	glm::dvec3 f(double t){ return o + r*glm::dvec3(cos(2 * glm::pi<float>()*t), sin(2 * glm::pi<float>()*t), 0.0); }
	glm::dvec3 df(double t){ return r*glm::dvec3(-sin(2 * glm::pi<float>()*t), cos(2 * glm::pi<float>()*t), 0.0); }
	glm::dvec3 ddf(double t){ return -r*glm::dvec3(cos(2 * glm::pi<float>()*t), sin(2 * glm::pi<float>()*t), 0.0); }
	glm::dvec3 dddf(double t){ return -r*glm::dvec3(-sin(2 * glm::pi<float>()*t), cos(2 * glm::pi<float>()*t), 0.0); }

private:
	glm::dvec3 o;
	double r;
};

class Spiral : public Curve
{
public:
	Spiral(glm::dvec3 p, double x) :o(p), a(x){}

	glm::dvec3 dnf(double t, unsigned n){ switch (n){ case 1: return df(t); case 2: return ddf(t); case 3: return dddf(t); default: throw UnsupportedDerivativeOrder(n); } }

	glm::dvec3 f(double t){ return o + a*t*glm::dvec3(cos(t), sin(t), 1); }
	glm::dvec3 df(double t){ return glm::dvec3(a*cos(t) - a*t*sin(t), a*sin(t) + a*t*cos(t), a); }
	glm::dvec3 ddf(double t){ return glm::dvec3(-a*sin(t) - a*t*cos(t) - a*sin(t), a*cos(t) - a*t*sin(t) + a*cos(t), 0); }
	glm::dvec3 dddf(double t){ return glm::dvec3(-a*cos(t) + a*t*sin(t) - a*cos(t) - a*cos(t), -a*sin(t) - a*t*cos(t) - a*sin(t) - a*sin(t), 0); }
private:
	glm::dvec3 o;
	double a;
};
