#include "MathFunctions.hpp"
#include "GeomInvariant.hpp"
#include <iostream> //temp

allBinomial<MAX_BINOMIAL_ORDER> binomialCoeff;


const float lobatto_roots[10] = {0, 0.040233045899999986, 0.13061306745, 0.2610375251, 0.41736052115, 0.58263947885, 0.7389624749, 0.86938693255, 0.9597669541, 1};
const float lobatto_coeffs[10] = {0.0111111111, 0.0666529954, 0.112444671, 0.1460213418, 0.1637698806, 0.1637698806, 0.1460213418, 0.112444671, 0.0666529954, 0.0111111111};

double integrate(std::function<double(Curve&,double)> func, Curve::Ptr c, double a, double b)
{
	double s = 0.0;

	for(int i = 0; i < 10; ++i)
	{
    double root = (b-a)*lobatto_roots[i]+a;
		double val = func(*c, root);
		s += val * val * GeomInv::v(*c, root) * GeomInv::v(*c, root) * lobatto_coeffs[i] * (b-a);
    //std::cerr << "x=" << root << std::endl << "dK=" << val << std::endl << "v=" << GeomInv::v(*c, root) << std::endl;
	}

	return sqrt(s);
}
