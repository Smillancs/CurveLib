#include "GeomInvariant.hpp"

#include "Utils.hpp"

namespace GeomInv{
	
	double frenet_coordinate(Curve& c, double t, unsigned xyz, unsigned order);

	glm::dvec3 e(Curve& c, double t)
	{ 
		auto normed = (c.dnf(t,1)!=glm::dvec3(0,0,0))?
			glm::normalize(c.dnf(t,1)):
			glm::dvec3(0,0,0);
		return normed;
	}
	
	glm::dvec3 b(Curve& c, double t, unsigned limit)
	{
		unsigned i = 2;
		while( glm::cross( c.dnf(t, 1), c.dnf(t, i) ) == glm::dvec3(0,0,0) && i < limit ) ++i;
		auto crossed = glm::cross(c.dnf(t, 1), c.dnf(t, i));
		auto normed = (crossed != glm::dvec3(0,0,0))?
			glm::normalize(crossed):
			glm::dvec3(0,0,0);
		return normed;
	}

	glm::dvec3 n(Curve& c, double t){ return glm::cross(b(c, t), e(c, t)); }

	double K(Curve& c, double t)
	{
		double x1 = frenet_coordinate(c, t, 1, 1);
		double y2 = frenet_coordinate(c, t, 2, 2);
		return y2 / (x1 * x1);
	}
	double T(Curve& c, double t)
	{
		double x1 = frenet_coordinate(c, t, 1, 1);
		double y2 = frenet_coordinate(c, t, 2, 2);
		double z3 = frenet_coordinate(c, t, 3, 3);
		return (y2 == 0) ? 0.0 : z3 / (x1 * y2);
	}


	double frenet_coordinate(Curve& c, double t, unsigned xyz, unsigned order)
	{
		switch(xyz)
		{
			case 1: return glm::dot(c.dnf(t,order), e(c,t));
			case 2: return glm::dot(c.dnf(t,order), n(c,t));
			case 3: return glm::dot(c.dnf(t,order), b(c,t));
			default: throw(Exception("The Frenet frame has only 3 dimensions, "+std::to_string(xyz)+" cannot be accepted."));
		}
	}
}
