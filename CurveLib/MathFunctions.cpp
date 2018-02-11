#include "MathFunctions.hpp"
#include "GeomInvariant.hpp"
#include <iostream> //temp

allBinomial<MAX_BINOMIAL_ORDER> binomialCoeff;


// Lobatto-quadrature
const double lobatto4_roots[4] = {0, 0.2763932023, 0.7236067977, 1};
const double lobatto4_coeffs[4] = {0.0833333333, 0.41666666665, 0.41666666665, 0.0833333333};
const double lobatto8_roots[8] = {0, 0.06412992575, 0.2041499093, 0.39535039105, 0.60464960895, 0.7958500907, 0.93587007425, 1};
const double lobatto8_coeffs[8] = {0.01785714285, 0.10535211355, 0.1705613462, 0.2062293973, 0.2062293973, 0.1705613462, 0.10535211355, 0.01785714285};
const double lobatto16_roots[16] = {0, 0.0152159769, 0.0503997335, 0.1039958541, 0.1738056486, 0.2569702891, 0.35008476555, 0.44933686325, 0.55066313675, 0.64991523445, 0.7430297109, 0.8261943514, 0.8960041459, 0.9496002665, 0.9847840231, 1};
const double lobatto16_coeffs[16] = {0.00416666665, 0.0290149465, 0.04469684865, 0.06212769105, 0.0770134904, 0.08874595665, 0.0968450119, 0.10097915405, 0.10097915405, 0.0968450119, 0.08874595665, 0.0770134904, 0.06212769105, 0.04469684865, 0.0290149465, 0.00416666665};
const double lobatto32_roots[32] = {0, 0.00369553305, 0.0123526548, 0.0258575808, 0.0440750305, 0.066823762, 0.09387763415, 0.12496775305, 0.1597851222, 0.19798370645, 0.2391838686, 0.2829761414, 0.32892529675, 0.3765746706, 0.4254507016, 0.4750676375, 0.5249323625, 0.5745492984, 0.6234253294, 0.67107470325, 0.7170238586, 0.7608161314, 0.80201629355, 0.8402148778, 0.87503224695, 0.90612236585, 0.933176238, 0.9559249695, 0.9741424192, 0.9876473452, 0.99630446695, 1};
const double lobatto32_coeffs[32] = {0.0010080645, 0.00619905325, 0.0110997764, 0.0158875677, 0.02051710075, 0.02494263565, 0.0291202486, 0.0330084386, 0.0365685698, 0.0397652628, 0.04256674895, 0.04494518645, 0.04687693775, 0.04834280445, 0.04932821825, 0.04982338575, 0.04982338575, 0.04932821825, 0.04834280445, 0.04687693775, 0.04494518645, 0.04256674895, 0.0397652628, 0.0365685698, 0.0330084386, 0.0291202486, 0.02494263565, 0.02051710075, 0.0158875677, 0.0110997764, 0.00619905325, 0.0010080645};
const double lobatto10_roots[10] = {0, 0.0402330459, 0.13061306745, 0.2610375251, 0.41736052115, 0.58263947885, 0.7389624749, 0.86938693255, 0.9597669541, 1};
const double lobatto10_coeffs[10] = {0.0111111111, 0.0666529954, 0.112444671, 0.1460213418, 0.1637698806, 0.1637698806, 0.1460213418, 0.112444671, 0.0666529954, 0.0111111111};

double lobatto_root(unsigned deg, unsigned num)
{
  if(deg == 4) return lobatto4_roots[num];
  else if(deg == 8) return lobatto8_roots[num];
  else if(deg == 16) return lobatto16_roots[num];
  else if(deg == 32) return lobatto32_roots[num];
}
double lobatto_coeff(unsigned deg, unsigned num)
{
  if(deg == 4) return lobatto4_coeffs[num];
  else if(deg == 8) return lobatto8_coeffs[num];
  else if(deg == 16) return lobatto16_coeffs[num];
  else if(deg == 32) return lobatto32_coeffs[num];
}

double integral2(std::function<double(Curve&,double)> func, Curve::Ptr c, double a, double b, unsigned deg)
{
	double s = 0.0;

  for(int i = 0; i < deg; ++i)
  {
    double root = (b-a)*lobatto_root(deg,i)+a;
    double val = func(*c, root);
    double vel = GeomInv::v(*c, root);
    s += val * val * vel * vel * lobatto_coeff(deg,i) * (b-a);
    //std::cerr << "x=" << root << std::endl << "dK=" << val << std::endl << "v=" << GeomInv::v(*c, root) << std::endl;
  }
  std::cerr << deg << ' ' << sqrt(s) << std::endl;
	return sqrt(s);
}

double integral2_uniform(std::function<double(Curve&,double)> func, Curve::Ptr c, double a, double b, unsigned deg)
{
	double s = 0.0;

  for(int i = 0; i < deg; ++i)
  {
    double root = (b-a)*(i/static_cast<float>(deg))+a;
    double val = func(*c, root);
    double vel = GeomInv::v(*c, root);
    s += val * val * vel * vel * (b-a) / deg;
    //std::cerr << "x=" << root << std::endl << "dK=" << val << std::endl << "v=" << GeomInv::v(*c, root) << std::endl;
  }
  //std::cerr << deg << ' ' << sqrt(s) << std::endl;
	return sqrt(s);
}

double integrate(std::function<double(Curve&,double)> func, Curve::Ptr c, double a, double b)
{
  double s = 0.0;
  double sprev = -1.0;
  double eps = 1e-3;
  unsigned deg = 4;
  //std::cerr << c->about() << std::endl;
  while(deg <= 10000 && abs(s-sprev) > eps)
  {
    sprev = s;
    s = integral2_uniform(func, c, a, b, deg);
    deg *= 2;
  }
  return s;
}



double maxvalue(std::function<double(Curve&,double)> func, Curve::Ptr c, double a, double b, unsigned num)
{
  double mx = 0.f;

  for(int i = 1; i < num; ++i)
  {
    double val = abs(func(*c, i/static_cast<float>(num)));
    if(val > mx) mx = val;
  }
  //std::cerr << num << ' ' << mx << std::endl;

  return mx;
}

double infnorm(std::function<double(Curve&,double)> func, Curve::Ptr c, double a, double b)
{
  double s = 0.f;
  double sprev = -1.f;
  double eps = 1e-3;
  uint num = 100;
  while(num <= 10000 && abs(s-sprev) > eps)
  {
    sprev = s;
    s = maxvalue(func, c, a, b, num);
    num *= 10;
  }
  return s;
}
