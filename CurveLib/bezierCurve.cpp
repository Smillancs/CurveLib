﻿#include "bezierCurve.h"

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
	size_t n = cp.size()-1;
	for (size_t i = 0; i <= n; ++i)
		coeff.push_back(binomial(n, i));
	dirtyCoeff = false;
}

void BezierCurve::CheckCoeff()
{
	if (dirtyCoeff)
	{
		RegenerateCoeff();
		dirtyCoeff = false;
	}
}

void BezierCurve::addControlPoint(const Eigen::Vector3f _cp)
{
	cp.push_back(_cp);
	dirtyCoeff = true;
}

//Farin 52. 4.5 Derivates and the de Casteljau algorithm
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

//Farin 44-45. 4.1 Bernstein Polynomials
double BezierCurve::Bernstein(const size_t n, const size_t i, const double t) const
{
	return binomial(n,i)*pow(t,i)*pow(1.0-t,n-i);
}

//Farin 44-45. 4.1 Bernstein Polynomials
double BezierCurve::Bernstein(const size_t i, const double t) const
{
	auto n = cp.size()-1;
	return coeff[i]*pow(t, i)*pow(1.0-t, n-i);
}

//Farin 44-45. 4.1 Bernstein Polynomials
Eigen::Vector3f BezierCurve::BernsteinEval(const float t)
{
	CheckCoeff();
	Eigen::Vector3f temp = Eigen::Vector3f(0,0,0);
	size_t n = cp.size();
	for(size_t i = 0; i < n; ++i)
	{
		temp += Bernstein(i,t)*cp[i];
	}
	return mainCoeff*temp;
}

//Farin 50. (4.4 Higher order derivates)
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

//Farin 49->51. 4.4 Higher order derivates
Eigen::Vector3f BezierCurve::ForwardDiff(const size_t order, const size_t idx)
{
	Eigen::Vector3f delta = {0,0,0};
	int alternate = (order%2)?-1:1;
	for (size_t j = 0; j <= order; ++j)
	{
		delta += binomial(order, j)*alternate*cp[idx+j];
		alternate*=-1;
	}
	return delta;
}

//Farin 49->51. 4.4 Higher order derivates
Eigen::Vector3f BezierCurve::Diff(const size_t order, const float t)
{
	CheckCoeff();
	Eigen::Vector3f eval = {0,0,0};
	size_t n = cp.size() - order -1;

	for (size_t j = 0; j <= n; ++j)
	{
		eval += ForwardDiff(order,j)*Bernstein(n,j,t);
	}
	eval*=factorial(cp.size()-1)/static_cast<double>(factorial(n));
	return eval;
}

//Farin 53 4.6 Subdivision
SubCurves BezierCurve::Subdivision(const double t)
{
	CheckCoeff();
	SubCurves dividedCmp;

	std::vector<Eigen::Vector3f> tempCp = cp;

	const size_t n = cp.size();

	for (size_t i = 0; i < n; ++i)
	{
		dividedCmp.first.cp[i] = Eigen::Vector3f(0, 0, 0);
		dividedCmp.second.cp[i] = Eigen::Vector3f(0, 0, 0);
		for (size_t j = 0; j <= i; ++j)
		{
			dividedCmp.first.cp[i] += cp[j]*Bernstein(i, j, t);
		}
		for (size_t j = 0; j <= n-i; ++j)
		{
			dividedCmp.second.cp[i] += cp[j]*Bernstein(n-i, j, t);
		}
	}

	return dividedCmp;
}

//Farin 64. 5.1 Degree elevation
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

//Farin 67->70. 5.4 Degree reduction
BezierCurve BezierCurve::Reduction()
{
	CheckCoeff();
	BezierCurve reduced;
	reduced.cp.clear();

	const size_t n = cp.size();
	std::vector<Eigen::Vector3f> biLR = {}, biRL = {};
	for (size_t i = 0; i < n-1; ++i)
	{
		biLR.push_back((n*cp[i]-i*biLR[i-1])/(n-i));
	}
	for (size_t i = n; i < 1; --i)
	{
		biRL.push_back((n*cp[i] - (n-i)*biRL[i+1])/i);
	}
	for (size_t i = 0; i < n-1; ++i)
	{
		double L = 0;
		for (size_t j = 0; j < i; ++j)
		{
			L += binomial(2*n,2*j);
		}
		L /= 1/(pow(2,2*n-1));
		reduced.addControlPoint((1-L)*biLR[i]+L*biRL[i]);
	}

	return reduced;
}


BezierCurve BezierCurve::ToParametric()
{
	return *this;
}

//Farin 71. 5.5 Nonparametric curves
BezierCurve::Explicit BezierCurve::ToExplicit(const double t)
{
	CheckCoeff();
	BezierCurve::Explicit eval;
	double n = cp.size();
	for (size_t i = 0; i < n; ++i)
	{
		double b = Bernstein(n, i, t);
		eval.first+=i/n*b;
		eval.second+=cp[i]*b;
	}
	eval.first/=static_cast<double>(n);
	return eval;
}

std::string BezierCurve::about()
{
	std::stringstream s;
	s << "Bezier curve with control points: ";
	Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "(", ")", "", "");
	for(auto it : cp)
		s << it.transpose().format(HeavyFmt) << "; ";
	return s.str();
}

#endif
