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

struct ControlPoint
{
	double x,y,z;
};

class BezierCurve : public Curve
{
public:
	BezierCurve();
	BezierCurve(const std::vector<Eigen::Vector3f> &_cp);
	BezierCurve(const std::vector<ControlPoint> &_cp);

	BezierCurve operator=(const BezierCurve & _other);

	std::vector<Eigen::Vector3f> deCasteljauEval(const float t, const size_t deg);
	Eigen::Vector3f BernsteinEval(const float t);
	
	BezierCurve Diff(const size_t order = 1, const bool returnAsHodo = false);

	SubCurves Subdivision(const double t);
	BezierCurve Elevation();
	BezierCurve Reduction();

	void ToExplicit();
	void ToParametric();

	void addControlPoint(const Eigen::Vector3f _cp);

	inline std::vector<Eigen::Vector3f> GetControlPoints() const { return cp; }


	virtual glm::dvec3 dnf(double t, unsigned n) override 
	{ 
		BezierCurve diff = Diff(n);
		auto diffEval = diff.BernsteinEval(t);

		return glm::dvec3(diffEval.x(),diffEval.y(),diffEval.z());
	}

	virtual glm::dvec3 f(double t) override
	{
		Eigen::Vector3f bez = BernsteinEval(t);
		return glm::dvec3(bez.x(), bez.y(), bez.z());
	}


private:
#ifdef GPGPU
	void initOGL();
	std::vector<GLuint> shaders;
	GLuint m_vaoID, m_vboID;
#endif

	std::vector<Eigen::Vector3f> cp;
	std::vector<double> coeff;
	bool dirtyCoeff;
	double mainCoeff;
	void RegenerateCoeff();
};


#endif //BEZIER_H
