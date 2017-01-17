#ifndef BEZIER_H
#define BEZIER_H
//#define GLEW_STATIC
#define GPGPU
#ifdef GPGPU
#include "GLUtils.h"
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform2.hpp>
#include <GL/glu.h>
#include <math.h>
#endif
#include <Eigen/Core>
#include "all_incl.h"

struct ControlPoint
{
	double x,y,z;
};

class BezierSquareSurf
{
public:
	BezierSquareSurf();
	BezierSquareSurf(const std::vector<std::vector<Eigen::Vector3f>> &_cp);
	BezierSquareSurf(const std::vector<std::vector<ControlPoint>> &_cp);

	BezierSquareSurf operator=(const BezierSquareSurf & _other);

	std::vector<Eigen::Vector3f> deCasteljauEval(const float t, const size_t deg);
	Eigen::Vector3f BernsteinEval(const float t);
	
	BezierSquareSurf Diff(const size_t order = 1, const bool returnAsHodo = false);

	//SubCurves Subdivision(const double t);
	BezierSquareSurf Elevation();
	BezierSquareSurf Reduction();

	void ToExplicit();
	void ToParametric();

	void addControlPoint(const Eigen::Vector3f _cp);

	inline std::vector<std::vector<Eigen::Vector3f>> GetControlPoints() const { return cp; }


	virtual glm::dvec3 dnf(double t, unsigned n)  
	{ 
		BezierSquareSurf diff = Diff(n);
		auto diffEval = diff.BernsteinEval(t);

		return glm::dvec3(diffEval.x(),diffEval.y(),diffEval.z());
	}

	virtual glm::dvec3 f(double t) 
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

	std::vector<std::vector<Eigen::Vector3f>> cp;
	std::vector<double> coeff;
	bool dirtyCoeff;
	double mainCoeff;
	void RegenerateCoeff();
};


#endif //BEZIER_H
