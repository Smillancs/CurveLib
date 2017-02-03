#pragma once

#include "Curve.hpp"

namespace GeomInv
{

	glm::dvec3 e(Curve& c, double t);
	glm::dvec3 b(Curve& c, double t, unsigned limit = 5);
	glm::dvec3 n(Curve& c, double t);

	double K(Curve& c, double t);
	double dK(Curve& c, double t);
	double ddK(Curve& c, double t);
	
	double T(Curve& c, double t);
	double dT(Curve& c, double t);

}
