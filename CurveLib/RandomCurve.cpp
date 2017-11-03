#include "RandomCurve.hpp"

#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <ctime>
#include "MathFunctions.hpp"


std::vector<double> randomCoeffList(int n)
{
  std::vector<double> v(n);
  for(auto& it : v) it = rand()/(RAND_MAX/2.0)-1.0;
  return v;
}

constexpr double transformationCoeff(int n, int j, int k)
{
  double m = 0;
  for(int i=std::max(0,j+k-n); i<=std::min(j,k); ++i)
  {
    m += ((k+i)&1?-1:1)*binomialCoeff(k,i)*binomialCoeff(k,i)*binomialCoeff(n-k,j-i);
  }
  m /= binomialCoeff(n,j);
  return m;
}

std::vector<double> LegendreToBernstein(std::vector<double> l, int deg)
{
  std::vector<double> b(deg+1);
  for(int j = 0; j < l.size(); ++j)
    for(int k = 0; k < l.size(); ++k)
    {
      b[j] += transformationCoeff(deg, j, k)*l[k];
    }
  return b;
}

Curve::Ptr RandomCurve(int deg, int dim)
{
  srand(time(0));
  std::vector<glm::vec3> controlPoints(deg+1);
  std::vector<double> x = LegendreToBernstein(randomCoeffList(deg+1), deg);
  std::vector<double> y = LegendreToBernstein(randomCoeffList(deg+1), deg);
  std::vector<double> z = LegendreToBernstein(randomCoeffList(deg+1), deg);
  int i = 0;
  for(auto& it : controlPoints)
  {
    it.x = x[i];
    it.y = y[i];
    it.z = dim>2 ? z[i] : 0;
    ++i;
  }
  Curve::Ptr c = Curve::Ptr(new BezierCurve(controlPoints));
  return c;
}
