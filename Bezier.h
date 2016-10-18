#ifndef BEZIER_H
#define BEZIER_H
#include <Eigen\Core>
#include "all_incl.h"

class BezierCurve
{
public:
	BezierCurve();
	
	BezierCurve(const std::vector<Eigen::Vector2f> &_cp);

	BezierCurve operator=(const BezierCurve & _other);

	std::vector<Eigen::Vector2f> deCasteljauEval(const float t, const size_t deg);
	Eigen::Vector2f BernsteinEval(const float t);
	
	BezierCurve Diff(const size_t order = 1);
	double DiffEval(const double t);

	BezierCurve Subdivision();
	BezierCurve Elevation();
	BezierCurve Reduction();

	void ToExplicit();
	void ToParametric();

	void addControlPoint(const Eigen::Vector2f _cp);

	inline std::vector<Eigen::Vector2f> GetControlPoints() const { return cp; }
private:
	std::vector<Eigen::Vector2f> cp;
	std::vector<double> coeff;
	bool dirtyCoeff;
	double mainCoeff;
	void RegenerateCoeff();
};


#endif //BEZIER_H