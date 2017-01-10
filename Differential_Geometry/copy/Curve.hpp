#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>

#include <math.h>

class Curve
{
public:
	virtual ~Curve(){}

	virtual glm::dvec3 f(double t) = 0;

	virtual glm::dvec3 dnf(double t, unsigned n) = 0;

	// These are intended to be removed later, if not needed.
	virtual glm::dvec3 df(double t){ return dnf(t, 1); }
	virtual glm::dvec3 ddf(double t){ return dnf(t, 2); }
	virtual glm::dvec3 dddf(double t){ return dnf(t, 3); }

};
