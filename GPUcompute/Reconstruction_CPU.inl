
#include "../CurveLib/MathFunctions.hpp"
#include "../CurveLib/Utils.hpp"

namespace std
{
  template<typename _Tp>
  struct array<_Tp,0>{
    _Tp& operator[](int){ throw("The impossible happened.");/*_Tp x; return x;*/}
  };
}

constexpr const std::array<float,4*4> pinverse1 = {1,0,0,0, 1,1/3.0,0,0, 0,0,1,-1/3.0, 0,0,1,0};
constexpr const std::array<float,6*6> pinverse2 = {1.0000,0.0000,0.0000,0,0,0,1.0000,0.2000,0.0000,0,0,0,1.0000,0.4000,0.0500,0,0,0,0,0,0,1.0000,-0.4000,0.0500,0,0,0,1.0000,-0.2000,-0.0000,0,0,0,1.0000,-0.0000,-0.0000};
constexpr const std::array<float,8*8> pinverse3 = {1.0000,-0.0000,-0.0000,0.0000,0,0,0,0,1.0000,0.1429,-0.0000,0.0000,0,0,0,0,1.0000,0.2857,0.0238,0.0000,0,0,0,0,1.0000,0.4286,0.0714,0.0048,0,0,0,0,0,0,0,0,1.0000,-0.4286,0.0714,-0.0048,0,0,0,0,1.0000,-0.2857,0.0238,-0.0000,0,0,0,0,1.0000,-0.1429,0.0000,-0.0000,0,0,0,0,1.0000,-0.0000,0.0000,-0.0000};

template<int continuity>
constexpr std::array<float,(2*continuity+2)*(2*continuity+2)> getPinverse()
{
  if constexpr(continuity == 1) return pinverse1;
  if constexpr(continuity == 2) return pinverse2;
  if constexpr(continuity == 3) return pinverse3;
}

template<int extra_points, size_t DoF>
std::array<glm::vec3,2> PointDerivatives1(const std::vector<typename Reconstruction<1,extra_points>::Input>& input, std::array<float,DoF> freedoms, int id, bool end = false)
{
  int idx = end ? 2*id+1 : 2*id;
  int fr = end ? 1 : 0;
  glm::vec3 d0 = freedoms[fr] * input[idx].e();
	return std::array<glm::vec3,2>{input[idx].p(), d0};
}
template<int extra_points, size_t DoF>
std::array<glm::vec3,3> PointDerivatives2(const std::vector<typename Reconstruction<2,extra_points>::Input>& input, std::array<float,DoF> freedoms, int id, bool end = false)
{
  int idx = end ? 2*id+1 : 2*id;
  int fr = end ? 1 : 0;
  int fr2 = end ? 3 : 2;
  glm::vec3 d0 = freedoms[fr] * input[idx].e();
  glm::vec3 dd0 = freedoms[fr2] * input[idx].e() + input[idx].K() * freedoms[fr] * freedoms[fr] * input[idx].n();
  return std::array<glm::vec3,3>{input[idx].p(), d0, dd0};
}
template<int extra_points, size_t DoF>
std::array<glm::vec3,4> PointDerivatives3(const std::vector<typename Reconstruction<3,extra_points>::Input>& input, std::array<float,DoF> freedoms, int id, bool end = false)
{
  int idx = end ? 2*id+1 : 2*id;
  int fr = end ? 1 : 0;
  int fr2 = end ? 3 : 2;
  int fr4 = end ? 5 : 4;
  glm::vec3 d0 = freedoms[fr] * input[idx].e();
  glm::vec3 dd0 = freedoms[fr2] * input[idx].e() + input[idx].K() * freedoms[fr] * freedoms[fr] * input[idx].n();
  glm::vec3 b0 = cross(input[idx].e(), input[idx].n());
  glm::vec3 ddd0 = freedoms[fr4] * input[idx].e() + (3*freedoms[fr]*freedoms[fr2]*input[idx].K() + freedoms[fr] * freedoms[fr] * freedoms[fr] * input[idx].dK(1)) * input[idx].n() + freedoms[fr] * freedoms[fr] * freedoms[fr] * input[idx].K() * input[idx].T() * b0;
  return std::array<glm::vec3,4>{input[idx].p(), d0, dd0, ddd0};
}

template<int continuity, int extra_points, size_t DoF>
std::array<glm::vec3,continuity+1> PointDerivatives(const std::vector<typename Reconstruction<continuity,extra_points>::Input>& input, std::array<float,DoF> freedoms, int id, bool end = false)
{
  if constexpr(continuity == 1) return PointDerivatives1<extra_points,DoF>(input, freedoms, id, end);
  if constexpr(continuity == 2) return PointDerivatives2<extra_points,DoF>(input, freedoms, id, end);
  if constexpr(continuity == 3) return PointDerivatives3<extra_points,DoF>(input, freedoms, id, end);
}

template<int continuity, int extra_points, size_t DoF, size_t point_num>
std::array<float, DoF> getFreedoms(Curve::Ptr c)
{
  std::array<float, DoF> freedoms;
  for(int i=0;i<continuity;++i)
  {
    freedoms[2*i] = glm::dot(c->dnf(0, i+1), GeomInv::e(*c, 0));
    freedoms[2*i+1] = glm::dot(c->dnf(1, i+1), GeomInv::e(*c, 1));
  }
  if constexpr(extra_points > 0)
  {
    std::vector<glm::vec3> v = std::dynamic_pointer_cast<BezierCurve>(c)->GetGlmControlPoints();
    std::array<glm::vec3,point_num> controlPoints;
    std::copy(v.begin(), v.end(), controlPoints.begin());
    for(int i=0;i<extra_points;++i)
    {
      freedoms[2*continuity+3*i]   = controlPoints[continuity+1+i].x;
      freedoms[2*continuity+3*i+1] = controlPoints[continuity+1+i].y;
      freedoms[2*continuity+3*i+2] = controlPoints[continuity+1+i].z;
    }
  }
  return freedoms;
}

template<int continuity, int extra_points, size_t point_num>
std::array<glm::vec3,point_num> calculateControlPoints(std::array<glm::vec3,continuity+1> start, std::array<glm::vec3,continuity+1> end, std::array<glm::vec3,extra_points> extras)
{
  std::array<float,(2*continuity+2)*3> pointdata;
  for(int i=0;i < continuity+1;++i)
  {
    pointdata[3*i+0] = start[i].x;
    pointdata[3*i+1] = start[i].y;
    pointdata[3*i+2] = start[i].z;
  }
  for(int i=continuity+1;i < 2*continuity+2;++i)
  {
    pointdata[3*i+0] = end[i-continuity-1].x;
    pointdata[3*i+1] = end[i-continuity-1].y;
    pointdata[3*i+2] = end[i-continuity-1].z;
  }
  std::array<float,(2*continuity+2)*3> controlpoints = matmul<continuity>(getPinverse<continuity>(), pointdata);
  std::array<glm::vec3,point_num> controlPoints;
  for(int i=0;i < continuity+1;++i)
    controlPoints[i] = glm::vec3(controlpoints[3*i],controlpoints[3*i+1],controlpoints[3*i+2]);
  for(int i=0;i < extra_points;++i)
    controlPoints[continuity+1+i] = extras[i];
  for(int i=0;i < continuity+1;++i)
    controlPoints[continuity+1+extra_points+i] = glm::vec3(controlpoints[3*(continuity+1)+3*i],controlpoints[3*(continuity+1)+3*i+1],controlpoints[3*(continuity+1)+3*i+2]);
  return controlPoints;
}

template<int continuity, int extra_points, size_t DoF, size_t point_num>
std::array<float,DoF> optimizeNewton(std::function<double(Curve&,double)> func, const std::vector<typename Reconstruction<continuity,extra_points>::Input>& input, std::array<float,DoF> start, int id, bool opt_points, bool infnorm, bool& no_move)
{
  const int ITERATIONS = 10;
  float eps = 1e-3;
  float min_dist = length(input[2*id].p() - input[2*id+1].p()) / 100.0;
  std::array<glm::vec3,point_num> points;
  std::array<float,DoF> freedoms;

  std::array<float,2*continuity> d;
  std::array<float,2*continuity*2*continuity> dd;
  std::array<float,2*continuity> newton_step;
  std::array<float,3*extra_points> _d;
  std::array<float,3*extra_points*3*extra_points> _dd;
  std::array<float,3*extra_points> _newton_step;
  std::array<std::array<float,DoF>,ITERATIONS+1> _positions;
  std::array<typename Reconstruction<continuity,extra_points>::Result_cpu,ITERATIONS+1> positions;

  for(int i=0;i < DoF;++i) freedoms[i] = start[i];

	for(int i=0;i < ITERATIONS;++i)
	{
    if(!opt_points){
      std::array<glm::vec3,extra_points> extras;
      for(int j=0;j < extra_points;++j)
        extras[j] = glm::vec3(freedoms[2*continuity+3*j+0],freedoms[2*continuity+3*j+1],freedoms[2*continuity+3*j+2]);
      std::array<float,9> integrals;
      for(int x=0;x < 2*continuity;++x) for(int y=0;y <= x;++y){
        for(int z=0;z < 9;++z){
          if((x == y && (z != 3 && z != 4 && z!= 5)) || (z == 1 || z == 7 )) continue;
          std::array<float,DoF> freedoms_local;
          for(int j=0;j < DoF;++j)
          {
            if( j == x ) freedoms_local[j] = freedoms[j]+(int(z)%3-1)*eps;
            else if( j == y ) freedoms_local[j] = freedoms[j]+(int(z)/3-1)*eps;
            else freedoms_local[j] = freedoms[j];
          }
          points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives<continuity,extra_points,DoF>(input, freedoms_local, id, false), PointDerivatives<continuity,extra_points,DoF>(input, freedoms_local, id, true), extras);
          /*if(infnorm)
      		  integrals[z] = infnorm(func, Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end()))));
          else*/
      		  integrals[z] = integrate(func, Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end()))));

          if(x == 0 && y == 0 && z == 4)
          {
            for(int j=0;j < point_num;++j)
      			   positions[i].first[j] = points[j];
      			positions[i].second = integrals[4];
            _positions[i] = freedoms;
            //std::cerr << BezierCurve(std::vector<glm::vec3>(points.begin(), points.end())).about() << std::endl;
          }

        }
        d[x] = (integrals[5]-integrals[3])/(2*eps);
        if(x == y) dd[x*2*continuity+y] = (integrals[5] - 2*integrals[4] + integrals[3])/(eps*eps);
        else{
          dd[x*2*continuity+y] = ((integrals[8]-integrals[6])-(integrals[2]-integrals[0]))/(4*eps*eps);
          dd[y*2*continuity+x] = ((integrals[8]-integrals[6])-(integrals[2]-integrals[0]))/(4*eps*eps);
        }
      }
      newton_step = linsolve(dd,d);
      //std::cerr << "Norm: " << integrals[4] << std::endl;
      /*printArray(d);
      printArray(dd);
      printArray(newton_step);*/
      for(int x=0;x < 2*continuity;++x)
        freedoms[x] -= newton_step[x];

      /*if(freedoms[0] < min_dist) freedoms[0] = min_dist;
      if(freedoms[1] < min_dist) freedoms[1] = min_dist;*/

      float n = 0.0;
      for(int j=0;j < 2*continuity;++j)
      {
        n += d[j]*d[j];
      }
      if(n < eps) break;
    }
    else{
      std::array<float,9> integrals;
      for(int x=2*continuity;x < DoF;++x) for(int y=2*continuity;y <= x;++y){
        for(int z=0;z < 9;++z){
          if((x == y && (z != 3 && z != 4 && z!= 5)) || (z == 1 || z == 7 )) continue;
          std::array<float,DoF> freedoms_local;
          for(int j=0;j < DoF;++j)
          {
            if( j == x ) freedoms_local[j] = freedoms[j]+(int(z)%3-1)*eps;
            else if( j == y ) freedoms_local[j] = freedoms[j]+(int(z)/3-1)*eps;
            else freedoms_local[j] = freedoms[j];
          }
          std::array<glm::vec3,extra_points> extras;
          for(int j=0;j < extra_points;++j)
            extras[j] = glm::vec3(freedoms_local[2*continuity+3*j+0],freedoms_local[2*continuity+3*j+1],freedoms_local[2*continuity+3*j+2]);
          points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives<continuity,extra_points,DoF>(input, freedoms_local, id, false), PointDerivatives<continuity,extra_points,DoF>(input, freedoms_local, id, true), extras);
          /*if(infnorm)
      		  integrals[z] = infnorm(func, Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end()))));
          else*/
      		  integrals[z] = integrate(func, Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end()))));

          if(x == 2*continuity && y == 2*continuity && z == 4)
          {
            for(int j=0;j < point_num;++j)
      			   positions[i].first[j] = points[j];
      			positions[i].second = integrals[4];
            _positions[i] = freedoms;
            std::cerr << BezierCurve(std::vector<glm::vec3>(points.begin(), points.end())).about() << std::endl;
          }
        }
        _d[-2*continuity+x] = (integrals[5]-integrals[3])/(2*eps);
        if(x == y) _dd[(-2*continuity+x)*3*extra_points+(-2*continuity+y)] = (integrals[5] - 2*integrals[4] + integrals[3])/(eps*eps);
        else{
          _dd[(-2*continuity+x)*3*extra_points+(-2*continuity+y)] = ((integrals[8]-integrals[6])-(integrals[2]-integrals[0]))/(eps*eps);
          _dd[(-2*continuity+y)*3*extra_points+(-2*continuity+x)] = ((integrals[8]-integrals[6])-(integrals[2]-integrals[0]))/(eps*eps);
        }
      }
      _newton_step = linsolve(_dd,_d);
      if constexpr(extra_points > 0)
      {
        std::cerr << "Norm: " << integrals[4] << std::endl;
        printArray(freedoms);
        printArray(_d);
        printArray(_dd);
        printArray(_newton_step);
      }
      for(int x=0;x < 3*extra_points;++x)
        freedoms[x+2*continuity] -= _newton_step[x];

      /*if(freedoms[0] < min_dist) freedoms[0] = min_dist;
      if(freedoms[1] < min_dist) freedoms[1] = min_dist;*/

      float n = 0.0;
      for(int j=0;j < 3*extra_points;++j)
      {
        n += _d[j]*_d[j];
      }
      if(n < eps) break;
    }
  }

  std::array<glm::vec3,extra_points> extras;
  for(int j=0;j < extra_points;++j)
    extras[j] = glm::vec3(freedoms[2*continuity+3*j+0],freedoms[2*continuity+3*j+1],freedoms[2*continuity+3*j+2]);
  points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives<continuity,extra_points,DoF>(input, freedoms, id, false), PointDerivatives<continuity,extra_points,DoF>(input, freedoms, id, true), extras);
  for(int j=0;j < point_num;++j)
     positions[ITERATIONS].first[j] = points[j];
  /*if(infnorm)
  	positions[ITERATIONS].norm = infnorm(func, Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end()))));
  else*/
    positions[ITERATIONS].second = integrate(func, Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end()))));
  _positions[ITERATIONS] = freedoms;

  int minimum = 0;
  for(int i=1;i<=ITERATIONS;++i)
  {
  	if(!isnan(positions[i].second) && positions[i].second < positions[minimum].second) minimum = i;
  }
  std::array<float,DoF> _current;
  _current = _positions[minimum];
  //current = positions[minimum];
  float n = 0.f;
  for(int i=0;i < DoF;++i)
  {
    n += pow(start[i]-_current[i], 2);
  }
  n = sqrt(n);
  no_move = n < eps;
  return _current;
}


template <int continuity, int extra_points>
std::vector<typename Reconstruction<continuity,extra_points>::Result_cpu> Reconstruction<continuity,extra_points>::optimize_cpu(std::function<double(Curve&,double)> func, const std::vector<typename Reconstruction<continuity,extra_points>::Input>& input, const std::shared_ptr<std::vector<float>>& debugInfo)
{
  static const size_t DoF = 2*continuity+3*extra_points;
  static const size_t point_num = 2*continuity+2+extra_points;
  std::array<glm::vec3,point_num> points;
  std::array<float,DoF> _start;
  std::array<float,DoF> _current;
  typename Reconstruction<continuity,extra_points>::Result current;
  std::vector<typename Reconstruction<continuity,extra_points>::Result_cpu> output(input.size()/2);

  for(int id = 0; id < input.size()/2; ++id)
  {
    for(int i=0;i < 2*continuity;++i) _start[i] = glm::length(input[2*id].p() - input[2*id+1].p()) / 5;
    for(int i=2*continuity;i < DoF;++i) _start[i] = 0;
    std::array<glm::vec3,extra_points> extras;
    for(int j=0;j < extra_points;++j)
      extras[j] = glm::vec3(_start[2*continuity+3*j+0],_start[2*continuity+3*j+1],_start[2*continuity+3*j+2]);
    points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives<continuity,extra_points,DoF>(input, _start, id, false), PointDerivatives<continuity,extra_points,DoF>(input, _start, id, true), extras);
    glm::vec3 point_a = points[continuity];
    glm::vec3 point_b = points[point_num-continuity-1];
    for(int i=0; i < extra_points; ++i)
    {
      glm::vec3 p = (static_cast<float>(extra_points-i)*point_a + static_cast<float>(i+1)*point_b) / static_cast<float>(extra_points+1);
      _start[2*continuity+3*i+0] = p.x;
      _start[2*continuity+3*i+1] = p.y;
      _start[2*continuity+3*i+2] = p.z;
    }
    _current = _start;

    bool finished = false;
    int c = 0;
      _current = optimizeLBFGS<DoF, continuity, extra_points>(func, input, _start, id);
    /*if(_current[0] < glm::length(input[2*id].p() - input[2*id+1].p()) / 100.0) _current[0] = glm::length(input[2*id].p() - input[2*id+1].p()) / 100.0;
    if(_current[1] < glm::length(input[2*id].p() - input[2*id+1].p()) / 100.0) _current[1] = glm::length(input[2*id].p() - input[2*id+1].p()) / 100.0;*/
    std::cerr << "Done" << std::endl;
    for(int j=0;j < extra_points;++j)
      extras[j] = glm::vec3(_current[2*continuity+3*j+0],_current[2*continuity+3*j+1],_current[2*continuity+3*j+2]);
    points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives<continuity,extra_points,DoF>(input, _current, id, false), PointDerivatives<continuity,extra_points,DoF>(input, _current, id, true), extras);
    //if(id==1) for(int i=0;i < 12;++i) dump.data[i] = positions[min].points[i/3][i%3];

    for(int j=0;j < point_num;++j)
      output[id].first[j] = points[j];
    output[id].second = integrate(func, Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end()))));
  }
  return output;
}

template <int continuity, int extra_points>
std::vector<typename Reconstruction<continuity,extra_points>::Result_cpu> Reconstruction<continuity,extra_points>::optimize_cpu_alt(std::function<double(Curve&,double)> func, const std::vector<typename Reconstruction<continuity,extra_points>::Input>& input, const std::shared_ptr<std::vector<float>>& debugInfo)
{
  static const size_t DoF = 2*continuity+3*extra_points;
  static const size_t point_num = 2*continuity+2+extra_points;
  std::array<glm::vec3,point_num> points;
  std::array<float,2*continuity> _start1;
  std::array<float,DoF> _start;
  std::array<float,2*continuity> _current1;
  std::array<float,DoF> _current;
  typename Reconstruction<continuity,extra_points>::Result current;
  std::vector<typename Reconstruction<continuity,extra_points>::Result_cpu> output(input.size()/2);

  for(int id = 0; id < input.size()/2; ++id)
  {
    for(int i=0;i < 2*continuity;++i) _start1[i] = glm::length(input[2*id].p() - input[2*id+1].p()) / 5;
    bool finished = false;

    _current1 = optimizeNewton<continuity,0,DoF-3*extra_points,point_num-extra_points>(func, input, _start1, id, false, false, finished);


    if constexpr(extra_points > 0)
    {
      auto points = calculateControlPoints<continuity,0,point_num-extra_points>(PointDerivatives<continuity,0,DoF-3*extra_points>(input, _current1, id, false), PointDerivatives<continuity,0,DoF-3*extra_points>(input, _current1, id, true), std::array<glm::vec3,0>());

      Curve::Ptr c = Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end())));

      BezierCurve bc;
      bc = *std::dynamic_pointer_cast<BezierCurve>(c);
      for(size_t i = 0; i < extra_points; ++i)
      {
        bc = bc.Elevation();
        //std::cerr << "Degree elevation: " << std::endl;
        //std::cerr << bc.about() << std::endl;
        //std::cerr << integrate(func, c) << std::endl;
      }
      c = Curve::Ptr(new BezierCurve(bc));

      _start = getFreedoms<continuity, extra_points, 2*continuity+3*extra_points, 2*continuity+2+extra_points>(c);

      _current = optimizeNewton<continuity,extra_points,DoF,point_num>(func, input, _start, id, true, false, finished);

    }
    else
      _current = _current1;

    /*if(_current[0] < glm::length(input[2*id].p() - input[2*id+1].p()) / 100.0) _current[0] = glm::length(input[2*id].p() - input[2*id+1].p()) / 100.0;
    if(_current[1] < glm::length(input[2*id].p() - input[2*id+1].p()) / 100.0) _current[1] = glm::length(input[2*id].p() - input[2*id+1].p()) / 100.0;*/
    //std::cerr << "Done" << std::endl;
    std::array<glm::vec3,extra_points> extras;
    for(int j=0;j < extra_points;++j)
      extras[j] = glm::vec3(_current[2*continuity+3*j+0],_current[2*continuity+3*j+1],_current[2*continuity+3*j+2]);
    points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives<continuity,extra_points,DoF>(input, _current, id, false), PointDerivatives<continuity,extra_points,DoF>(input, _current, id, true), extras);
    //if(id==1) for(int i=0;i < 12;++i) dump.data[i] = positions[min].points[i/3][i%3];

    for(int j=0;j < point_num;++j)
      output[id].first[j] = points[j];
    output[id].second = integrate(func, Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end()))));
  }
  return output;
}

glm::vec3 operator*(double a, glm::vec3 v){ return static_cast<float>(a)*v; }


// unused code with dlib
/*#include "../dlib/optimization.h"
#include "../dlib/global_optimization.h"


typedef dlib::matrix<double,0,1> column_vector;

template<int extra_points, size_t DoF>
std::array<glm::vec3,2> PointDerivatives1_bfgs(const std::vector<typename Reconstruction<1,extra_points>::Input>& input, column_vector freedoms, int id, bool end = false)
{
  int idx = end ? 2*id+1 : 2*id;
  int fr = end ? 1 : 0;
  glm::vec3 d0 = freedoms(fr) * input[idx].e();
	return std::array<glm::vec3,2>{input[idx].p(), d0};
}
template<int extra_points, size_t DoF>
std::array<glm::vec3,3> PointDerivatives2_bfgs(const std::vector<typename Reconstruction<2,extra_points>::Input>& input, column_vector freedoms, int id, bool end = false)
{
  int idx = end ? 2*id+1 : 2*id;
  int fr = end ? 1 : 0;
  int fr2 = end ? 3 : 2;
  glm::vec3 d0 = freedoms(fr) * input[idx].e();
  glm::vec3 dd0 = freedoms(fr2) * input[idx].e() + input[idx].K() * freedoms(fr) * freedoms(fr) * input[idx].n();
  return std::array<glm::vec3,3>{input[idx].p(), d0, dd0};
}
template<int extra_points, size_t DoF>
std::array<glm::vec3,4> PointDerivatives3_bfgs(const std::vector<typename Reconstruction<3,extra_points>::Input>& input, column_vector freedoms, int id, bool end = false)
{
  int idx = end ? 2*id+1 : 2*id;
  int fr = end ? 1 : 0;
  int fr2 = end ? 3 : 2;
  int fr4 = end ? 5 : 4;
  glm::vec3 d0 = freedoms(fr) * input[idx].e();
  glm::vec3 dd0 = freedoms(fr2) * input[idx].e() + input[idx].K() * freedoms(fr) * freedoms(fr) * input[idx].n();
  glm::vec3 b0 = cross(input[idx].e(), input[idx].n());
  glm::vec3 ddd0 = freedoms(fr4) * input[idx].e() + (3*freedoms(fr)*freedoms(fr2)*input[idx].K() + freedoms(fr) * freedoms(fr) * freedoms(fr) * input[idx].dK(1)) * input[idx].n() + freedoms(fr) * freedoms(fr) * freedoms(fr) * input[idx].K() * input[idx].T() * b0;
  return std::array<glm::vec3,4>{input[idx].p(), d0, dd0, ddd0};
}

template<int continuity, int extra_points, size_t DoF>
std::array<glm::vec3,continuity+1> PointDerivatives_bfgs(const std::vector<typename Reconstruction<continuity,extra_points>::Input>& input, column_vector freedoms, int id, bool end = false)
{
  if constexpr(continuity == 1) return PointDerivatives1_bfgs<extra_points,DoF>(input, freedoms, id, end);
  if constexpr(continuity == 2) return PointDerivatives2_bfgs<extra_points,DoF>(input, freedoms, id, end);
  if constexpr(continuity == 3) return PointDerivatives3_bfgs<extra_points,DoF>(input, freedoms, id, end);
}


template<int continuity, int extra_points>
struct functor{
  const typename std::vector<typename Reconstruction<continuity,extra_points>::Input> input;
  int id;
  std::function<double(Curve&,double)> func;

  static const int point_num = 2*continuity+2+extra_points;
  static const int DoF = 2*continuity+3*extra_points;

  functor(const typename std::vector<typename Reconstruction<continuity,extra_points>::Input> in, int i, std::function<double(Curve&,double)> f)
    : input(in), id(i), func(f) {}

  double operator() (const column_vector& m) const
  {
    std::array<glm::vec3,extra_points> extras;
    if constexpr(extra_points > 0)
    {
      for(int j=0;j < extra_points;++j)
        extras[j] = glm::vec3(m(2*continuity+3*j+0),m(2*continuity+3*j+1),m(2*continuity+3*j+2));
    }
    auto points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives_bfgs<continuity,0,DoF-3*extra_points>(input, m, id, false), PointDerivatives_bfgs<continuity,0,DoF-3*extra_points>(input, m, id, true), extras);
    Curve::Ptr c = Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end())));
    return integrate(func, Curve::Ptr(c));
  }
};

template <int continuity, int extra_points>
std::vector<typename Reconstruction<continuity,extra_points>::Result_cpu> Reconstruction<continuity,extra_points>::optimize_cpu_bfgs(std::function<double(Curve&,double)> func, const std::vector<typename Reconstruction<continuity,extra_points>::Input>& input, const std::shared_ptr<std::vector<float>>& debugInfo)
{
  static const size_t DoF = 2*continuity+3*extra_points;
  static const size_t point_num = 2*continuity+2+extra_points;
  std::array<glm::vec3,point_num> points;
  column_vector _start(DoF);
  column_vector _current(DoF);
  typename Reconstruction<continuity,extra_points>::Result current;
  std::vector<typename Reconstruction<continuity,extra_points>::Result_cpu> output(input.size()/2);

  for(int id = 0; id < input.size()/2; ++id)
  {
    for(int i=0;i < 2*continuity;++i) _start(i) = glm::length(input[2*id].p() - input[2*id+1].p()) / 5;
    for(int i=2*continuity;i < DoF;++i) _start(i) = 0;
    std::array<glm::vec3,extra_points> extras;
    for(int j=0;j < extra_points;++j)
      extras[j] = glm::vec3(0);
    points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives_bfgs<continuity,extra_points,DoF>(input, _start, id, false), PointDerivatives_bfgs<continuity,extra_points,DoF>(input, _start, id, true), extras);
    glm::vec3 point_a = points[continuity];
    glm::vec3 point_b = points[point_num-continuity-1];
    for(int i=0; i < extra_points; ++i)
    {
      glm::vec3 p = (static_cast<float>(extra_points-i)*point_a + static_cast<float>(i+1)*point_b) / static_cast<float>(extra_points+1);
      _start(2*continuity+3*i+0) = p.x;
      _start(2*continuity+3*i+1) = p.y;
      _start(2*continuity+3*i+2) = p.z;
    }
    _current = _start;

    bool finished = false;
    int c = 0;

    functor<continuity, extra_points> f(input, id, func);
    
    dlib::find_min_using_approximate_derivatives(dlib::bfgs_search_strategy(),
                                           dlib::objective_delta_stop_strategy(1e-20),
                                           f, _current, -1);

    std::cerr << "Done" << std::endl;
    for(int j=0;j < extra_points;++j)
      extras[j] = glm::vec3(_current(2*continuity+3*j+0),_current(2*continuity+3*j+1),_current(2*continuity+3*j+2));
    points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives_bfgs<continuity,extra_points,DoF>(input, _current, id, false), PointDerivatives_bfgs<continuity,extra_points,DoF>(input, _current, id, true), extras);
    //if(id==1) for(int i=0;i < 12;++i) dump.data[i] = positions[min].points[i/3][i%3];

    for(int j=0;j < point_num;++j)
      output[id].first[j] = points[j];
    output[id].second = integrate(func, Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end()))));
  }
  return output;
}


template<int continuity, int extra_points>
struct functor_frenet_only{
  const typename std::vector<typename Reconstruction<continuity,extra_points>::Input> input;
  int id;
  std::function<double(Curve&,double)> func;
  std::array<glm::vec3, extra_points> extras;

  static const int point_num = 2*continuity+2+extra_points;
  static const int DoF = 2*continuity+3*extra_points;

  functor_frenet_only(const typename std::vector<typename Reconstruction<continuity,extra_points>::Input> in, int i, std::function<double(Curve&,double)> f, std::array<glm::vec3, extra_points> p)
    : input(in), id(i), func(f), extras(p) {}

  double operator() (const column_vector& m) const
  {
    auto points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives_bfgs<continuity,0,DoF-3*extra_points>(input, m, id, false), PointDerivatives_bfgs<continuity,0,DoF-3*extra_points>(input, m, id, true), extras);
    Curve::Ptr c = Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end())));
    return integrate(func, Curve::Ptr(c));
  }
};

template<int continuity, int extra_points>
struct functor_points_only{
  const typename std::vector<typename Reconstruction<continuity,extra_points>::Input> input;
  int id;
  std::function<double(Curve&,double)> func;
  column_vector frenet;

  static const int point_num = 2*continuity+2+extra_points;
  static const int DoF = 2*continuity+3*extra_points;

  functor_points_only(const typename std::vector<typename Reconstruction<continuity,extra_points>::Input> in, int i, std::function<double(Curve&,double)> f, column_vector x)
    : input(in), id(i), func(f), frenet(x) {}

  double operator() (const column_vector& m) const
  {
    std::array<glm::vec3,extra_points> extras;
    if constexpr(extra_points > 0)
    {
      for(int j=0;j < extra_points;++j)
        extras[j] = glm::vec3(m(3*j+0),m(3*j+1),m(3*j+2));
    }
    auto points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives_bfgs<continuity,0,DoF-3*extra_points>(input, frenet, id, false), PointDerivatives_bfgs<continuity,0,DoF-3*extra_points>(input, frenet, id, true), extras);
    Curve::Ptr c = Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end())));
    return integrate(func, Curve::Ptr(c));
  }
};

template <int continuity, int extra_points>
std::vector<typename Reconstruction<continuity,extra_points>::Result_cpu> Reconstruction<continuity,extra_points>::optimize_cpu_bfgs_elev(std::function<double(Curve&,double)> func, const std::vector<typename Reconstruction<continuity,extra_points>::Input>& input, const std::shared_ptr<std::vector<float>>& debugInfo)
{
  static const size_t DoF = 2*continuity+3*extra_points;
  static const size_t point_num = 2*continuity+2+extra_points;
  std::array<glm::vec3,point_num> points;
  std::array<glm::vec3, extra_points> extras;
  column_vector _start1(2*continuity, 1);
  column_vector _current1(2*continuity, 1);
  typename Reconstruction<continuity,extra_points>::Result current;
  std::vector<typename Reconstruction<continuity,extra_points>::Result_cpu> output(input.size()/2);

  for(int id = 0; id < input.size()/2; ++id)
  {
    for(int i=0;i < 2*continuity;++i) _start1(i) = glm::length(input[2*id].p() - input[2*id+1].p()) / 5;
    _current1 = _start1;

    bool finished = false;

    functor_frenet_only<continuity, 0> f1(input, id, func, std::array<glm::vec3,0>());
    
    dlib::find_min_using_approximate_derivatives(dlib::bfgs_search_strategy(),
                                           dlib::objective_delta_stop_strategy(1e-7),
                                           f1, _current1, -1);

    std::array<glm::vec3, point_num - extra_points> points1 = calculateControlPoints<continuity,0,point_num-extra_points>(PointDerivatives_bfgs<continuity,extra_points,DoF>(input, _current1, id, false), PointDerivatives_bfgs<continuity,extra_points,DoF>(input, _current1, id, true), std::array<glm::vec3, 0>());
    BezierCurve bc (std::vector<glm::vec3>(points1.begin(), points1.end()));
    std::cerr << bc.about() << std::endl;
    std::cerr << integrate(func, Curve::Ptr(new BezierCurve(bc))) << std::endl;

    if constexpr(extra_points > 0)
    {
      column_vector _start2(3*extra_points, 1);
      column_vector _current2(3*extra_points, 1);
      for(int i=0;i<extra_points;++i)
      {
        bc = bc.Elevation();
      }
      Curve::Ptr c = Curve::Ptr(new BezierCurve(bc));
      std::array<float, DoF> arr = getFreedoms<continuity, extra_points, DoF, point_num>(c);

      for(int i=0;i<DoF;++i)
      {
        if(i < 2*continuity)
          _start1(i) = arr[i];
        else
          _start2(i-2*continuity) = arr[i];
      }

      functor_points_only<continuity, extra_points> f2(input, id, func, _start1);
      _current1 = _start1;
      _current2 = _start2;

      std::cerr << _start1 << std::endl;

      dlib::find_min_using_approximate_derivatives(dlib::bfgs_search_strategy(),
                                            dlib::objective_delta_stop_strategy(1e-7),
                                            f2, _current2, -1);

      std::cerr << f2(_start2) << std::endl;
      std::cerr << f2(_current2) << std::endl;
      for(int j=0;j < extra_points;++j)
        extras[j] = glm::vec3(_current2(3*j+0),_current2(3*j+1),_current2(3*j+2));
    }
    std::cerr << "Done" << std::endl;
    points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives_bfgs<continuity,extra_points,DoF>(input, _current1, id, false), PointDerivatives_bfgs<continuity,extra_points,DoF>(input, _current1, id, true), extras);
    //if(id==1) for(int i=0;i < 12;++i) dump.data[i] = positions[min].points[i/3][i%3];

    for(int j=0;j < point_num;++j)
      output[id].first[j] = points[j];
    output[id].second = integrate(func, Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end()))));
  }
  return output;
}*/


#include "../alglib/optimization.h"


typedef alglib::real_1d_array column_vec;

template<int extra_points, size_t DoF>
std::array<glm::vec3,2> PointDerivatives1_alglib(const std::vector<typename Reconstruction<1,extra_points>::Input>& input, column_vec freedoms, int id, bool end = false)
{
  int idx = end ? 2*id+1 : 2*id;
  int fr = end ? 1 : 0;
  glm::vec3 d0 = freedoms(fr) * input[idx].e();
	return std::array<glm::vec3,2>{input[idx].p(), d0};
}
template<int extra_points, size_t DoF>
std::array<glm::vec3,3> PointDerivatives2_alglib(const std::vector<typename Reconstruction<2,extra_points>::Input>& input, column_vec freedoms, int id, bool end = false)
{
  int idx = end ? 2*id+1 : 2*id;
  int fr = end ? 1 : 0;
  int fr2 = end ? 3 : 2;
  glm::vec3 d0 = freedoms(fr) * input[idx].e();
  glm::vec3 dd0 = freedoms(fr2) * input[idx].e() + input[idx].K() * freedoms(fr) * freedoms(fr) * input[idx].n();
  return std::array<glm::vec3,3>{input[idx].p(), d0, dd0};
}
template<int extra_points, size_t DoF>
std::array<glm::vec3,4> PointDerivatives3_alglib(const std::vector<typename Reconstruction<3,extra_points>::Input>& input, column_vec freedoms, int id, bool end = false)
{
  int idx = end ? 2*id+1 : 2*id;
  int fr = end ? 1 : 0;
  int fr2 = end ? 3 : 2;
  int fr4 = end ? 5 : 4;
  glm::vec3 d0 = freedoms(fr) * input[idx].e();
  glm::vec3 dd0 = freedoms(fr2) * input[idx].e() + input[idx].K() * freedoms(fr) * freedoms(fr) * input[idx].n();
  glm::vec3 b0 = cross(input[idx].e(), input[idx].n());
  glm::vec3 ddd0 = freedoms(fr4) * input[idx].e() + (3*freedoms(fr)*freedoms(fr2)*input[idx].K() + freedoms(fr) * freedoms(fr) * freedoms(fr) * input[idx].dK(1)) * input[idx].n() + freedoms(fr) * freedoms(fr) * freedoms(fr) * input[idx].K() * input[idx].T() * b0;
  return std::array<glm::vec3,4>{input[idx].p(), d0, dd0, ddd0};
}

template<int continuity, int extra_points, size_t DoF>
std::array<glm::vec3,continuity+1> PointDerivatives_alglib(const std::vector<typename Reconstruction<continuity,extra_points>::Input>& input, column_vec freedoms, int id, bool end = false)
{
  if constexpr(continuity == 1) return PointDerivatives1_alglib<extra_points,DoF>(input, freedoms, id, end);
  if constexpr(continuity == 2) return PointDerivatives2_alglib<extra_points,DoF>(input, freedoms, id, end);
  if constexpr(continuity == 3) return PointDerivatives3_alglib<extra_points,DoF>(input, freedoms, id, end);
}


template<int continuity, int extra_points>
struct functor_alglib{
  const typename std::vector<typename Reconstruction<continuity,extra_points>::Input> input;
  int id;
  std::function<double(Curve&,double)> func;

  static const int point_num = 2*continuity+2+extra_points;
  static const int DoF = 2*continuity+3*extra_points;

  functor_alglib(const typename std::vector<typename Reconstruction<continuity,extra_points>::Input> in, int i, std::function<double(Curve&,double)> f)
    : input(in), id(i), func(f) {}

  void operator() (const column_vec& m, double& res, void* ptr)
  {
    std::array<glm::vec3,extra_points> extras;
    if constexpr(extra_points > 0)
    {
      for(int j=0;j < extra_points;++j)
        extras[j] = glm::vec3(m(2*continuity+3*j+0),m(2*continuity+3*j+1),m(2*continuity+3*j+2));
    }
    auto points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives_alglib<continuity,0,DoF-3*extra_points>(input, m, id, false), PointDerivatives_alglib<continuity,0,DoF-3*extra_points>(input, m, id, true), extras);
    Curve::Ptr c = Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end())));
    res = integrate(func, Curve::Ptr(c));
  }
};

template<int continuity, int extra_points>
void f_ptr_alglib(const column_vec& m, double& res, void* ptr)
{
  functor_alglib<continuity,extra_points>* func = (functor_alglib<continuity,extra_points>*)ptr;
  func->operator()(m, res, ptr);
}

template <int continuity, int extra_points>
std::vector<typename Reconstruction<continuity,extra_points>::Result_cpu> Reconstruction<continuity,extra_points>::optimize_cpu_alglib(std::function<double(Curve&,double)> func, const std::vector<typename Reconstruction<continuity,extra_points>::Input>& input, const std::shared_ptr<std::vector<float>>& debugInfo)
{
  static const size_t DoF = 2*continuity+3*extra_points;
  static const size_t point_num = 2*continuity+2+extra_points;
  std::array<glm::vec3,point_num> points;
  std::array<double, DoF> arr;
  column_vec _start;
  _start.setcontent(DoF,&arr[0]);
  column_vec _current;
  _current.setcontent(DoF,&arr[0]);
  typename Reconstruction<continuity,extra_points>::Result current;
  std::vector<typename Reconstruction<continuity,extra_points>::Result_cpu> output(input.size()/2);

  for(int id = 0; id < input.size()/2; ++id)
  {
    for(int i=0;i < 2*continuity;++i) _start(i) = glm::length(input[2*id].p() - input[2*id+1].p()) / 5;
    for(int i=2*continuity;i < DoF;++i) _start(i) = 0;
    std::array<glm::vec3,extra_points> extras;
    for(int j=0;j < extra_points;++j)
      extras[j] = glm::vec3(0);
    points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives_alglib<continuity,extra_points,DoF>(input, _start, id, false), PointDerivatives_alglib<continuity,extra_points,DoF>(input, _start, id, true), extras);
    glm::vec3 point_a = points[continuity];
    glm::vec3 point_b = points[point_num-continuity-1];
    for(int i=0; i < extra_points; ++i)
    {
      glm::vec3 p = (static_cast<float>(extra_points-i)*point_a + static_cast<float>(i+1)*point_b) / static_cast<float>(extra_points+1);
      _start(2*continuity+3*i+0) = p.x;
      _start(2*continuity+3*i+1) = p.y;
      _start(2*continuity+3*i+2) = p.z;
    }
    _current = _start;

    bool finished = false;
    int c = 0;

    functor_alglib<continuity, extra_points> f(input, id, func);
    //real_1d_array x = "[0,0]";
    double epsg = 0.0000000001;
    double epsf = 0;
    double epsx = 0;
    double diffstep = 1.0e-6;
    alglib::ae_int_t maxits = 0;
    alglib::minlbfgsstate state;
    alglib::minlbfgsreport rep;

    alglib::minlbfgscreatef(1, _start, diffstep, state);
    alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxits);
    alglib::minlbfgsoptimize(state, f_ptr_alglib<continuity,extra_points>, 0, (void*) &f);
    alglib::minlbfgsresults(state, _current, rep);

    if(_current[0] < 0.1) _current[0] = 0.1;
    if(_current[1] < 0.1) _current[1] = 0.1;

    //printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
    //printf("%s\n", _current.tostring(2).c_str()); // EXPECTED: [-3,3]

    /*if(_current[0] < glm::length(input[2*id].p() - input[2*id+1].p()) / 100.0) _current[0] = glm::length(input[2*id].p() - input[2*id+1].p()) / 100.0;
    if(_current[1] < glm::length(input[2*id].p() - input[2*id+1].p()) / 100.0) _current[1] = glm::length(input[2*id].p() - input[2*id+1].p()) / 100.0;*/
    std::cerr << "Done" << std::endl;
    for(int j=0;j < extra_points;++j)
      extras[j] = glm::vec3(_current(2*continuity+3*j+0),_current(2*continuity+3*j+1),_current(2*continuity+3*j+2));
    points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives_alglib<continuity,extra_points,DoF>(input, _current, id, false), PointDerivatives_alglib<continuity,extra_points,DoF>(input, _current, id, true), extras);
    //if(id==1) for(int i=0;i < 12;++i) dump.data[i] = positions[min].points[i/3][i%3];

    for(int j=0;j < point_num;++j)
      output[id].first[j] = points[j];
    output[id].second = integrate(func, Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end()))));
  }
  return output;
}


/*#include "NLPInterfacePack_ExampleNLPBanded.hpp"
#include "MoochoPack_MoochoSolver.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_GlobalMPISession.hpp"

template <int continuity, int extra_points>
std::vector<typename Reconstruction<continuity,extra_points>::Result_cpu> Reconstruction<continuity,extra_points>::optimize_cpu_moocho(std::function<double(Curve&,double)> func, const std::vector<typename Reconstruction<continuity,extra_points>::Input>& input, const std::shared_ptr<std::vector<float>>& debugInfo)
{
  static const size_t DoF = 2*continuity+3*extra_points;
  static const size_t point_num = 2*continuity+2+extra_points;


  namespace mp  = MoochoPack;
  namespace nlpip = NLPInterfacePack;
  using mp::MoochoSolver;
  using nlpip::NLP;
  using nlpip::ExampleNLPBanded;
  typedef nlpip::size_type  size_type;
  typedef nlpip::value_type value_type;
  using Teuchos::CommandLineProcessor;
  bool success = true;

  //Teuchos::GlobalMPISession mpiSession;

  try {

    MoochoSolver solver;
  
    //
    // Get options from the command line
    //
    
    bool     show_options = false;
    int      nD = 1;
    int      nI = 1;
    int      bw = 1;
    int      mI = 0;
    double   xo = 1;
    bool     nlp_selects_basis = true;
    double   xDl = -NLP::infinite_bound();
    double   xDu = +NLP::infinite_bound(); 
    double   xIl = -NLP::infinite_bound(); 
    double   xIu = +NLP::infinite_bound();
    /*double   xDl = -1000.;
    double   xDu = +1000.;
    double   xIl = -1000.;
    double   xIu = +1000.;* /

    int      mU = 0;
    double   hl = -NLP::infinite_bound();
    double   hu = +NLP::infinite_bound();
    double   diag_scal = 1.0;
    double   diag_vary = 1.0;
    bool     sym_basis = false;
    double   f_offset  = 0.0;
    double   co        = 0.0;
    bool     ignore_constraints = false;
    
    //
    // Create and solve the NLP
    //

    ExampleNLPBanded
      nlp(nD,nI,bw,mU,mI,xo,xDl,xDu,xIl,xIu,hl,hu
        ,nlp_selects_basis,diag_scal,diag_vary
        ,sym_basis,f_offset,co,ignore_constraints
        );

    solver.set_nlp( Teuchos::rcp(&nlp,false) );

    const MoochoSolver::ESolutionStatus
      solution_status = solver.solve_nlp();
    
    //return solution_status;

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cout,success)

  //return MoochoSolver::SOLVE_RETURN_EXCEPTION;
}*/



template<size_t DoF, int continuity, int extra_points>
float optimizationTarget(std::function<double(Curve&,double)> func, const std::vector<typename Reconstruction<continuity,extra_points>::Input>& input, std::array<float,DoF> vals, uint id)
{
  const int point_num = 2*continuity+2+extra_points;
  std::array<glm::vec3,point_num> points = calculateControlPoints<continuity,extra_points,point_num>(PointDerivatives<continuity,extra_points,DoF>(input,vals,id,false), PointDerivatives<continuity,extra_points,DoF>(input,vals,id,true), std::array<glm::vec3,0>{});
  return integrate(func, Curve::Ptr(new BezierCurve(std::vector<glm::vec3>(points.begin(), points.end()))));
}

template<size_t DoF>
float optimizationTarget_dummy(std::array<float,DoF> vals)
{
  return pow(vals[0]-2,2) + pow(vals[1]-3,2) - 1;
}

template<size_t DoF>
float optimizationTarget_rosen(std::array<float,DoF> vals)
{
  return pow(1 - vals[0], 2) + 100 * pow(vals[1] - pow(vals[0], 2), 2);
}

#define ITERATIONS 5

template<size_t DoF>
std::array<float,DoF> get(std::array<float,ITERATIONS*DoF> t, int k)
{
  if(k < 0)
  {
    std::array<float,DoF> v;
    for(int i = 0; i < DoF; ++i)
    {
      v[i] = 0;
    }
    return v;
  }
  std::array<float,DoF> v;
  for(int i = 0; i < DoF; ++i)
  {
    v[i] = t[k*DoF+i];
  }
  return v;
}

template<size_t DoF>
std::array<float,DoF> add(std::array<float,DoF> t, std::array<float,DoF> a)
{
  for(int i = 0; i < DoF; ++i)
  {
    t[i] += a[i];
  }
  return t;
}

template<size_t DoF>
std::array<float,DoF> mult(std::array<float,DoF> t, float a)
{
  for(int i = 0; i < DoF; ++i)
  {
    t[i] *= a;
  }
  return t;
}

template<size_t DoF>
float dot(std::array<float,DoF> v1, std::array<float,DoF> v2)
{
  float s = 0;
  for(int i = 0; i < DoF; ++i)
  {
    s += v1[i] * v2[i];
  }
  return s;
}

template<size_t DoF>
std::array<float,DoF*DoF> diad(std::array<float,DoF> v1, std::array<float,DoF> v2)
{
  std::array<float,DoF*DoF> m;
  for(int i = 0; i < DoF; ++i) for(int j = 0; j < DoF; ++j)
  {
    m[i*DoF+j] = v1[i] * v2[j];
  }
  return m;
}

template<size_t DoF>
std::array<float,DoF> matmul(std::array<float,DoF*DoF> mat, std::array<float,DoF> v)
{
  std::array<float,DoF> res{0};
  for(int i=0;i < DoF;++i)
    for(int k=0;k < DoF;++k)
      res[i] += mat[i*DoF+k]*v[k];
  return res;
}


template<size_t DoF, int continuity, int extra_points>
std::array<float,DoF> optimizeLBFGS(std::function<double(Curve&,double)> func, const std::vector<typename Reconstruction<continuity,extra_points>::Input>& input, std::array<float,DoF> start, uint id)
{
  const float lower_bound = 0.1;
  const float eps = 1e-3;
  const int m = 5;

  std::array<float,ITERATIONS*DoF> x;
  std::array<float,ITERATIONS*DoF> g;
  std::array<float,ITERATIONS*DoF> s;
  std::array<float,ITERATIONS*DoF> y;
  std::array<float,ITERATIONS> ro;

  std::array<float,DoF> freedoms;
  for(int i=0;i < DoF;++i) freedoms[i] = start[i];
  for(int i=0;i < DoF;++i) x[i] = start[i];

  printArray(freedoms);

  for(volatile int i = 0; i < ITERATIONS; ++i)
  {
    std::array<float,DoF> d;
    float val = optimizationTarget<DoF,continuity,extra_points>(func, input, freedoms, id);
    for(int k = 0; k < DoF; ++k)
    {
      std::array<float,DoF> freedoms_fwd = freedoms;
      freedoms_fwd[k] += eps;
      float val_fwd = optimizationTarget<DoF,continuity,extra_points>(func, input, freedoms_fwd, id);
      std::array<float,DoF> freedoms_bck = freedoms;
      freedoms_bck[k] -= eps;
      float val_bck = optimizationTarget<DoF,continuity,extra_points>(func, input, freedoms_bck, id);

      d[k] = (val_fwd - val_bck) / (2 * eps);
    }
    for(int j = 0; j < DoF; ++j) g[i*DoF+j] = d[j];
    if(i > 0) for(int j = 0; j < DoF; ++j) y[(i-1)*DoF+j] = g[i*DoF+j] - g[(i-1)*DoF+j];
    if(i > 0) ro[i-1] = 1/dot(get<DoF>(y,i-1), get<DoF>(s,i-1));

    std::array<float,DoF> q;
    q = d;
    //printArray(d);

    std::array<float,DoF> alpha;
    for(int j = i-1; j >= std::max(0,i-m); --j)
    {
      float alpha_i = alpha[j-(i-m)] = ro[j]*dot(get<DoF>(s,j),q);
      q = add(q, mult(get<DoF>(y,j),-1*alpha_i));
    }
    //printArray(q);
    
    std::array<float,DoF*DoF> H;
    if(i > 0)
    {
      float h = 1 / dot(get<DoF>(y,i-1), get<DoF>(y,i-1));
      H = diad(get<DoF>(y,i-1), mult(get<DoF>(s,i-1), h));
    }
    else
    {
      for(int a=0;a<DoF;++a) for(int b=0;b<DoF;++b)
        if(a==b) H[a*DoF+b] = 1;
        else H[a*DoF+b] = 0;
    }

    std::array<float,DoF> z;
    z = matmul(H, q);

    for(int j = std::max(0,i-m); j <= i-1; ++j)
    {
      float beta_i = ro[j]*dot(get<DoF>(y,j),z);
      z = add(z,mult(get<DoF>(s,j), alpha[j-(i-m)] - beta_i));
    }

    for(int k = 0; k < DoF; ++k) freedoms[k] -= z[k];
    printArray(freedoms);

    for(int j = 0; j < DoF; ++j) x[(i+1)*DoF+j] = freedoms[j];
    for(int j = 0; j < DoF; ++j) s[i*DoF+j] = x[(i+1)*DoF+j] - x[i*DoF+j];
    std::cerr << i << std::endl;
  }
  return get<DoF>(x, ITERATIONS-1);
}

template <int continuity, int extra_points>
Curve::Ptr Reconstruction<continuity,extra_points>::createResultCurve(const typename Reconstruction<continuity,extra_points>::Result_cpu& res)
{
	std::vector<glm::vec3> cps(res.first.begin(), res.first.end());
	return Curve::Ptr(new BezierCurve(cps));
}
