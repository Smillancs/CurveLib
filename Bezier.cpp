#include "Bezier.h"

BezierCurve::BezierCurve() 
	: mainCoeff(1.0), cp({}), coeff({})
{
}

BezierCurve::BezierCurve(const std::vector<Eigen::Vector2f> &_cp)
	: cp(_cp), mainCoeff(1.0)
{
	RegenerateCoeff();
}


BezierCurve BezierCurve::operator=(const BezierCurve & _other)
{
	cp = _other.GetControlPoints();
	RegenerateCoeff();
	return *this;
}

void BezierCurve::RegenerateCoeff()
{
	coeff.clear();
	size_t n = cp.size();
	for (size_t i = 0; i < n; ++i)
		coeff.push_back(binomial(n, i));
	dirtyCoeff = false;
}

void BezierCurve::addControlPoint(const Eigen::Vector2f _cp)
{
	cp.push_back(_cp);
	dirtyCoeff = true;
}

std::vector<Eigen::Vector2f> BezierCurve::deCasteljauEval(const float t, const size_t deg)
{
	std::vector<Eigen::Vector2f> tempCp = cp;

	for (size_t j = 0; j < deg; ++j)
	{
		for (size_t i = 0; i < tempCp.size()-1; ++i)
		{
			tempCp[i] = (1.0f-t)*tempCp[i]+t*tempCp[i+1];
		}
		tempCp.pop_back();
	}

	return tempCp;
}

Eigen::Vector2f BezierCurve::BernsteinEval(const float t)
{
	Eigen::Vector2f temp = Eigen::Vector2f(0,0);
	size_t n = cp.size();
	for(size_t i = 0; i < n; ++i)
	{
		temp += coeff[i]*pow(t,i)*pow(1-t,n-i)*cp[i];
	}
	return mainCoeff*temp;
}

//todo cachelni az evalhoz
BezierCurve BezierCurve::Diff(const size_t order)
{
	BezierCurve hodo;
	size_t n = cp.size()-order;
	for(size_t i = 0; i < n; ++i)
	{
		Eigen::Vector2f lcp = {0,0};
		for (size_t j = 0; j < order; ++j)
		{
			lcp += binomial(order,j)*pow(-1,order-j)*cp[i+j];
		}
		hodo.addControlPoint(lcp);
	}
	hodo.RegenerateCoeff();
	hodo.mainCoeff = factorial(cp.size())/factorial(n);
	return hodo;
}

double BezierCurve::DiffEval(const double t)
{
	return 0.0;
}


BezierCurve BezierCurve::Subdivision()
{
	return BezierCurve();
}

BezierCurve BezierCurve::Elevation()
{
	return BezierCurve();
}

BezierCurve BezierCurve::Reduction()
{
	return BezierCurve();
}


void BezierCurve::ToExplicit()
{

}


void BezierCurve::ToParametric()
{

}
