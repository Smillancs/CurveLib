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

Eigen::Vector3f BezierCurve::BernsteinEval(const float t)
{
	Eigen::Vector3f temp = Eigen::Vector3f(0,0);
	size_t n = cp.size();
	for(size_t i = 0; i < n; ++i)
	{
		temp += coeff[i]*pow(t,i)*pow(1-t,n-i)*cp[i];
	}
	return mainCoeff*temp;
}

//todo cachelni az evalhoz
BezierCurve BezierCurve::Diff(const size_t order, const bool returnAsHodo)
{
	BezierCurve diffCoeff;
	size_t n = cp.size()-order;
	for(size_t i = 0; i < n; ++i)
	{
		Eigen::Vector3f lcp = {0,0};
		for (size_t j = 0; j < order; ++j)
		{
			lcp += binomial(order,j)*pow(-1,order-j)*cp[i+j];
		}
		diffCoeff.addControlPoint(lcp);
	}
	if (returnAsHodo)
	{
		auto tmpSub = diffCoeff.cp[0], tmpMain = diffCoeff.cp[1];
		for (size_t i = 1; i < diffCoeff.cp.size(); ++i)
		{
			auto tmp = diffCoeff.cp[i];
			diffCoeff.cp[i] = tmpMain - tmpSub;
			tmpSub = tmpMain, tmpMain = tmp;
		}
	}

	diffCoeff.RegenerateCoeff();
	diffCoeff.mainCoeff = factorial(cp.size())/factorial(n);

	return diffCoeff;
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