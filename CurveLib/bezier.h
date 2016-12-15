#ifndef BEZIER_H
#define BEZIER_H
#define GLEW_STATIC
//#define GPGPU
#ifdef GPGPU
#include "GLUtils.h"
#include <GL/glew.h>
#include <glm.hpp>
#include <gtc/matrix_transform.hpp>
#include <gtx/transform2.hpp>
#include <GL/GLU.h>
#include <math.h>
#endif
#include <Eigen\Core>
#include "all_incl.h"

class BezierCurve;

typedef std::pair<BezierCurve, BezierCurve> SubCurves;

class BezierCurve
{
public:
	BezierCurve();
	

	BezierCurve(const std::vector<Eigen::Vector2f> &_cp);

	BezierCurve operator=(const BezierCurve & _other);

	std::vector<Eigen::Vector2f> deCasteljauEval(const float t, const size_t deg);
	Eigen::Vector2f BernsteinEval(const float t);
	
	BezierCurve Diff(const size_t order = 1, const bool returnAsHodo = false);
	//double DiffEval(const double t);

	SubCurves Subdivision(const double t);
	BezierCurve Elevation();
	BezierCurve Reduction();

	void ToExplicit();
	void ToParametric();

	void addControlPoint(const Eigen::Vector2f _cp);

	inline std::vector<Eigen::Vector2f> GetControlPoints() const { return cp; }
private:
#ifdef GPGPU
	void initOGL();
	std::vector<GLuint> shaders;
	GLuint m_vaoID, m_vboID;
#endif

	std::vector<Eigen::Vector2f> cp;
	std::vector<double> coeff;
	bool dirtyCoeff;
	double mainCoeff;
	void RegenerateCoeff();
};


#endif //BEZIER_H