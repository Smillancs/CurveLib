#include "bezierCurve.h"

BezierCurve::BezierCurve()
	: mainCoeff(1.0)
{
}

BezierCurve::BezierCurve(const std::vector<Eigen::Vector3f> &_cp)
	: cp(_cp), mainCoeff(1.0)
{
	RegenerateCoeff();
}

BezierCurve::BezierCurve(const std::vector<glm::vec3> &_cp)
	: mainCoeff(1.0)
{
	for (auto &p : _cp)
	{
		cp.push_back(Eigen::Vector3f(p.x,p.y,p.z));
	}
	RegenerateCoeff();
}

#ifndef GPGPU


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

void BezierCurve::addControlPoint(const Eigen::Vector3f _cp)
{
	cp.push_back(_cp);
	dirtyCoeff = true;
}

std::vector<Eigen::Vector3f> BezierCurve::deCasteljauEval(const float t, const size_t deg)
{
	std::vector<Eigen::Vector3f> tempCp = cp;

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

double BezierCurve::Bernstein(const size_t n, const size_t i, const double t) const
{
	return binomial(n,i)*pow(t,i)*pow(1.0-t,n-i);
}

double BezierCurve::Bernstein(const size_t i, const double t) const
{
	auto n = cp.size();
	return coeff[i]*pow(t, i)*pow(1.0-t, n-i);
}

Eigen::Vector3f BezierCurve::BernsteinEval(const float t)
{
	Eigen::Vector3f temp = Eigen::Vector3f(0,0,0);
	size_t n = cp.size();
	for(size_t i = 0; i < n; ++i)
	{
		temp += Bernstein(i,t)*cp[i];
	}
	return mainCoeff*temp;
}

std::vector<Eigen::Vector3f> BezierCurve::Hodo(const size_t order)
{
	size_t n = cp.size()-order;
	std::vector<Eigen::Vector3f> coeff(n);
	for(size_t i = 0; i < n; ++i)
	{
		coeff.push_back(ForwardDiff(order,i));
	}
	auto tmpSub = coeff[0], tmpMain = coeff[1];
	for (size_t i = 1; i < coeff.size(); ++i)
	{
		auto tmp = coeff[i];
		coeff[i] = tmpMain - tmpSub;
		tmpSub = tmpMain, tmpMain = tmp;
	}

	return coeff;
}

Eigen::Vector3f BezierCurve::ForwardDiff(const size_t order, const size_t idx)
{
	Eigen::Vector3f delta = {0,0,0};
	int alternate = (order%2)?-1:1;
	for (size_t j = 0; j < order; ++j)
	{
		delta += binomial(order, j)*alternate*cp[idx+j];
		alternate*=-1;
	}
	return delta;
}


Eigen::Vector3f BezierCurve::Diff(const size_t order, const float t)
{
	Eigen::Vector3f eval;
	size_t n = cp.size() - order;

	for (size_t j = 0; j < n; ++j)
	{
		eval += ForwardDiff(order,j)*Bernstein(n,j,t);
	}

	return eval;
}

SubCurves BezierCurve::Subdivision(const double t)
{
	SubCurves dividedCmp;

	std::vector<Eigen::Vector3f> tempCp = cp;

	const size_t n = cp.size();

	for (size_t i = 0; i < n; ++i)
	{
		dividedCmp.first.cp[i] = Eigen::Vector3f(0, 0);
		dividedCmp.second.cp[i] = Eigen::Vector3f(0, 0);
		for (size_t j = 0; j <= i; ++j)
		{
			dividedCmp.first.cp[i] += cp[j]*binomial(i, j)*pow(t, j)*pow(1-t, n-j);
		}
		for (size_t j = 0; j <= n-i; ++j)
		{
			dividedCmp.second.cp[i] += cp[j]*binomial(n-i, j)*pow(t, j)*pow(1-t, n-j);
		}
	}

	return dividedCmp;
}

BezierCurve BezierCurve::Elevation()
{
	BezierCurve bez;
	bez.cp = {};
	size_t n = cp.size();
	cp.push_back(cp[0]);
	for (size_t i = 1; i < n; ++i)
	{
		bez.cp.push_back(i/(n+1.0)*cp[i-1] + (1.0 - i/(n+1.0))*cp[i]);
	}
	cp.push_back(cp[cp.size()-1]);
	return bez;
}

BezierCurve BezierCurve::Reduction()
{
	const size_t n = cp.size();
	std::vector<Eigen::Vector3f> biLR = {}, biRL = {};
	for (size_t i = 0; i < n-1; ++i)
	{
		biLR.push_back((n*cp[i]-i*biLR[i-1])/(n-1));
	}
	for (size_t i = n; i < 1; --i)
	{
		biRL.push_back((n*cp[i] - (n-i)*biRL[i+1])/i);
	}
	//todo lin komb csebisevvel farin 86
	return BezierCurve();
}


void BezierCurve::ToExplicit()
{
	//todo farin 87
}


void BezierCurve::ToParametric()
{

}//ez kb return this

#endif
