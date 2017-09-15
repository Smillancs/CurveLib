#include "Utils.hpp"

void exportCurveData(std::string filename, Curve::Ptr c, std::function<double(Curve&,double)> func, int n, int segmentNumber)
{
  std::ofstream out(filename.c_str(), std::ios::app);
  for(int i=0;i<=n;++i)
  {
    out << segmentNumber+(double)i/n << " " << func(*c,(double)i/n) << std::endl;
  }
  out.close();
}
