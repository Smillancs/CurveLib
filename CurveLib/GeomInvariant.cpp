#include "GeomInvariant.hpp"

#include "Utils.hpp"

namespace GeomInv{

	// For "C++ historical reasons" all these templates have to be listed here
	template <int order, int xyz>
	double frenet_coordinate(Curve& c, double t);

	template <>
	double frenet_coordinate<1,1>(Curve& c, double t)
	{
		return glm::dot(c.dnf(t,1), e(c,t));
	}

	template <>
	double frenet_coordinate<2,1>(Curve& c, double t)
	{
		return glm::dot(c.dnf(t,2), e(c,t));
	}

	template <>
	double frenet_coordinate<3,1>(Curve& c, double t)
	{
		return glm::dot(c.dnf(t,3), e(c,t));
	}

	template <>
	double frenet_coordinate<4,1>(Curve& c, double t)
	{
		return glm::dot(c.dnf(t,4), e(c,t));
	}

	template <>
	double frenet_coordinate<1,2>(Curve& c, double t)
	{
		return glm::dot(c.dnf(t,1), n(c,t));
	}

	template <>
	double frenet_coordinate<2,2>(Curve& c, double t)
	{
		return glm::dot(c.dnf(t,2), n(c,t));
	}

	template <>
	double frenet_coordinate<3,2>(Curve& c, double t)
	{
		return glm::dot(c.dnf(t,3), n(c,t));
	}

	template <>
	double frenet_coordinate<4,2>(Curve& c, double t)
	{
		return glm::dot(c.dnf(t,4), n(c,t));
	}

	template <>
	double frenet_coordinate<1,3>(Curve& c, double t)
	{
		return glm::dot(c.dnf(t,1), b(c,t));
	}

	template <>
	double frenet_coordinate<2,3>(Curve& c, double t)
	{
		return glm::dot(c.dnf(t,2), b(c,t));
	}

	template <>
	double frenet_coordinate<3,3>(Curve& c, double t)
	{
		return glm::dot(c.dnf(t,3), b(c,t));
	}

	template <>
	double frenet_coordinate<4,3>(Curve& c, double t)
	{
		return glm::dot(c.dnf(t,4), b(c,t));
	}

	template <int order, int xyz>
	double frenet_coordinate(Curve& c, double t)
	{
		throw(Exception("The Frenet frame has only 3 dimensions, "+std::to_string(xyz)+" cannot be accepted."));
	}

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

  double v(Curve& c, double t)
  {
    return glm::length(c.dnf(t,1));
  }

	double K(Curve& c, double t)
	{
		double x1 = frenet_coordinate<1,1>(c, t);
		double y2 = frenet_coordinate<2,2>(c, t);
		return y2 / (x1 * x1);
	}

	double dK(Curve& c, double t)
	{
		double x1 = frenet_coordinate<1,1>(c, t);
		double x2 = frenet_coordinate<2,1>(c, t);
		double y2 = frenet_coordinate<2,2>(c, t);
		double y3 = frenet_coordinate<3,2>(c, t);
		return (y3 * x1 - 3 * x2 * y2) / (x1 * x1 * x1 * x1);
	}

	double ddK(Curve& c, double t)
	{
		double x1 = frenet_coordinate<1,1>(c, t);
		double x2 = frenet_coordinate<2,1>(c, t);
		double x3 = frenet_coordinate<3,1>(c, t);
		double y2 = frenet_coordinate<2,2>(c, t);
		double y3 = frenet_coordinate<3,2>(c, t);
		double y4 = frenet_coordinate<4,2>(c, t);
		return (y4 * x1 * x1 - 5 * x2 * y3 * x1 + 12 * x2 * x2 * y2 - 4 * y2 * x3 * x1 - 3 * y2 * y2 * y2) / (x1 * x1 * x1 * x1 * x1 * x1);
	}

	double T(Curve& c, double t)
	{
		double x1 = frenet_coordinate<1,1>(c, t);
		double y2 = frenet_coordinate<2,2>(c, t);
		double z3 = frenet_coordinate<3,3>(c, t);
		return (y2 == 0) ? 0.0 : z3 / (x1 * y2);
	}

	double dT(Curve& c, double t)
	{
		double x1 = frenet_coordinate<1,1>(c, t);
		double y2 = frenet_coordinate<2,2>(c, t);
		double y3 = frenet_coordinate<3,2>(c, t);
		double z3 = frenet_coordinate<3,3>(c, t);
		double z4 = frenet_coordinate<4,3>(c, t);
		return (y2 == 0) ? 0.0 : (z4 * y2 - 2 * z3 * y3) / (x1 * x1 * y2 * y2);
	}

}
