#ifndef BEZIER_H
#define BEZIER_H
//#define GPGPU
#ifdef GPGPU
#include "GLUtils.h"
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform2.hpp>
#include <GL/GLU.h>
#include <math.h>
#endif
#include <Eigen/Core>
#include "all_incl.h"
#include "Curve.hpp"

class BezierCurve;

typedef std::pair<BezierCurve, BezierCurve> SubCurves;

class BezierCurve : public Curve
{
public:
	BezierCurve();
	BezierCurve(const std::vector<Eigen::Vector3f> &_cp);
	BezierCurve(const std::vector<glm::vec3> &_cp);

	BezierCurve operator=(const BezierCurve & _other);

	typedef std::pair<double, Eigen::Vector3f> Explicit;

	std::vector<Eigen::Vector3f> deCasteljauEval(const float t, const size_t deg);
	Eigen::Vector3f BernsteinEval(const float t);
	
	Eigen::Vector3f ForwardDiff(const size_t order, const size_t idx);
	Eigen::Vector3f Diff(const size_t order, const float t);
	std::vector<Eigen::Vector3f> Hodo(const size_t order);

	SubCurves Subdivision(const double t);
	BezierCurve Elevation();
	BezierCurve Reduction();

	Explicit ToExplicit(const double t);
	BezierCurve ToParametric();

	void addControlPoint(const Eigen::Vector3f _cp);

	std::vector<Eigen::Vector3f> GetControlPoints() const { return cp; }
	std::vector<glm::vec3> GetGlmControlPoints() const 
	{
		std::vector<glm::vec3> tempcp = {}; 
		for(auto &i : cp)
		{
			tempcp.push_back(glm::vec3(i.x(),i.y(),i.z()));
		}
		return tempcp; 
	}
	std::vector<double> GetBernsteinCoeffs() const { return coeff; }
	

	virtual glm::dvec3 dnf(double t, unsigned n) override 
	{ 
		Eigen::Vector3f diff = Diff(n,t);

		return glm::dvec3(diff.x(), diff.y(), diff.z());
	}

	virtual glm::dvec3 f(double t) override
	{
		Eigen::Vector3f bez = BernsteinEval(t);
		return glm::dvec3(bez.x(), bez.y(), bez.z());
	}
	
	virtual std::string about() override;


private:
#ifdef GPGPU
	void initOGL();
	std::vector<GLuint> shaders;
	GLuint m_vaoID, m_vboID;
#endif
	double Bernstein(const size_t n, const size_t i, const double t) const;
	double Bernstein(const size_t i, const double t) const;
	std::vector<Eigen::Vector3f> cp;
	std::vector<double> coeff;
	bool dirtyCoeff;
	double mainCoeff;
	void RegenerateCoeff();
	void CheckCoeff();

};


#endif //BEZIER_H
