
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
            std::cerr << BezierCurve(std::vector<glm::vec3>(points.begin(), points.end())).about() << std::endl;
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
      std::cerr << "Norm: " << integrals[4] << std::endl;
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
        /*printArray(_d);
        printArray(_dd);
        printArray(_newton_step);*/
      }
      for(int x=0;x < 3*extra_points;++x)
        freedoms[x+2*continuity] -= newton_step[x];

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
    if constexpr(extra_points > 0)
      _current = optimizeNewton<continuity,extra_points,DoF,point_num>(func, input, _start, id, true, false, finished);
    while(true)
    {
      _start = _current;

      _current = optimizeNewton<continuity,extra_points,DoF,point_num>(func, input, _start, id, false, false, finished);

      if(finished || c >= 2) break;

      if constexpr(extra_points > 0)
      {
         _start = _current;

         _current = optimizeNewton<continuity,extra_points,DoF,point_num>(func, input, _start, id, true, false, finished);

         if(finished || c >= 2) break;
      }
      c++;
    }
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
Curve::Ptr Reconstruction<continuity,extra_points>::createResultCurve(const typename Reconstruction<continuity,extra_points>::Result_cpu& res)
{
	std::vector<glm::vec3> cps(res.first.begin(), res.first.end());
	return Curve::Ptr(new BezierCurve(cps));
}
