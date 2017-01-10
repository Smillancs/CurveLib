#include "GeomInvariant.hpp"

namespace GeomInv{

	glm::dvec3 e(Curve& c, double t)
	{ 
		auto normed = (c.dnf(t,1)!=glm::dvec3(0,0,0))?
			glm::normalize(c.dnf(t,1)):
			glm::dvec3(0,0,0);
		return normed;
	}
	
	glm::dvec3 b(Curve& c, double t)
	{ 
		auto crossed = glm::cross(c.dnf(t, 1), c.dnf(t, 2));
		auto normed = (crossed!=glm::dvec3(0,0,0))?
			glm::normalize(crossed):
			glm::dvec3(0,0,0);
		return normed; 
	}

	glm::dvec3 n(Curve& c, double t){ return glm::cross(b(c, t), e(c, t)); }

	double K(Curve& c, double t){ return glm::dot(c.dnf(t,2), n(c, t)); }
	double T(Curve& c, double t){ return glm::dot(c.dnf(t,3), b(c, t)); }

}