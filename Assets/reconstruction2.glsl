 /* This file is generated from reconstruction.php with generate_reconstruction_glsl.sh
 * For permanent modifications, modify the php file */
#version 430

const int continuity = 2;
const int extra_points = 0;
const int point_num = (2*continuity+1)+1+extra_points;
const int max_length = point_num;
const int freedom = 4;

layout (local_size_x = 1, local_size_y = 1, local_size_z = 1) in;


struct ReconstructionData2
{
  vec3 p;
  vec3 e;
  vec3 n;
  float K;
};

struct status
{
  float freedoms[freedom];
};

struct Result
{
  vec3 points[point_num];
  float norm;
};

layout(std430, binding = 0) buffer inputBuffer
{
	ReconstructionData2 data[]; // will be paired i.e. startpoint, endpoint, startpoint, ...
} inBuf;

layout(std430, binding = 1) buffer destBuffer
{
  Result data[]; // data of one curve
} outBuf;

layout(std430, binding = 2) buffer debugInfo
{
	float data[100];
} dump;

#define _have_dump

#incl ../Assets/BezierEval.glsl

#incl ../Assets/GeomInvariant.glsl

subroutine float GeomInvariantEval(vec3 data[max_length], int point_num, float t);

subroutine ( GeomInvariantEval ) float curvature(vec3 data[max_length], int point_num, float t) {
    return K(data, point_num, t);
}
subroutine ( GeomInvariantEval ) float curvatureD(vec3 data[max_length], int point_num, float t) {
    return dK(data, point_num, t);
}
subroutine ( GeomInvariantEval ) float const_velocity(vec3 data[max_length], int point_num, float t) {
    return dot(Diff(data, 1, t, point_num), Diff(data, 2, t, point_num));
}

subroutine uniform GeomInvariantEval Eval;

uniform bool infnorm = false;

// Gauss–Legendre quadrature formula
const float legendre_coeffs[20] = float[20](0.00880700357, 0.0203007149, 0.03133602417, 0.04163837079, 0.05096505991, 0.05909726598, 0.06584431922, 0.07104805466, 0.07458649324, 0.07637669357, 0.07637669357, 0.07458649324, 0.07104805466, 0.06584431922, 0.05909726598, 0.05096505991, 0.04163837079, 0.03133602417, 0.0203007149, 0.00880700357);
const float legendre_roots[20]  = float[20](0.003435700407, 0.01801403636, 0.04388278587, 0.08044151409, 0.1268340468, 0.1819731596, 0.244566499, 0.3131469556, 0.3861070744, 0.4617367394, 0.5382632606, 0.6138929256, 0.6868530444, 0.755433501, 0.8180268404, 0.8731659532, 0.9195584859, 0.9561172141, 0.9819859636, 0.9965642996);
// Lobatto-quadrature
const float lobatto_roots[10] = float[10](0, 0.040233045899999986, 0.13061306745, 0.2610375251, 0.41736052115, 0.58263947885, 0.7389624749, 0.86938693255, 0.9597669541, 1);
const float lobatto_coeffs[10] = float[10](0.0111111111, 0.0666529954, 0.112444671, 0.1460213418, 0.1637698806, 0.1637698806, 0.1460213418, 0.112444671, 0.0666529954, 0.0111111111);

float integral2(vec3 points[max_length])
{
	float s = 0.f;

	for(int i = 0; i < 10; ++i)
	{
		float val = Eval(points, point_num, lobatto_roots[i]);
		s += val * val * v(points, point_num, lobatto_roots[i]) * v(points, point_num, lobatto_roots[i]) * lobatto_coeffs[i];
	}

	return sqrt(s);
}

float maxvalue(vec3 points[max_length])
{
  float mx = 0.f;

  for(int i = 1; i < 1000; ++i)
  {
    float val = abs(Eval(points, point_num, i/1000.f));
    if(val > mx) mx = val;
  }

  return mx;
}

float[(2*continuity+2)*3] matmul(float pinverse[(2*continuity+2)*(2*continuity+2)], float pointdata[(2*continuity+2)*3])
{
  float res[(2*continuity+2)*3];
  for(int i=0;i < (2*continuity+2);++i)
    for(int j=0;j < 3;++j)
    {
      res[3*i+j] = 0;
      for(int k=0;k < (2*continuity+2);++k)
        res[3*i+j] += pinverse[(2*continuity+2)*i+k]*pointdata[3*k+j];
    }
  return res;
}

float[4] linsolve(float mat[16], float vec[4])
{
  for(int k=0;k < 4;++k)
  {
    for(int i=k+1;i < 4;++i)
    {
      float f = mat[i*4+k] / mat[k*4+k];
      for(int j=k+1;j < 4;++j)
      {
        mat[i*4+j] -= mat[k*4+j] * f;
      }
      vec[i] -= vec[k]*f;
      mat[i*4+k] = 0;
    }
  }

  for(int i=4-1;i>=0;--i)
  {
    vec[i] /= mat[i*4+i];
    mat[i*4+i] = 1;
    for(int j=i-1;j>=0;--j)
    {
      vec[j] -= mat[j*4+i] * vec[i];
      mat[j*4+i] = 0;
    }
  }
  return vec;
}

const float pinverse[(2*continuity+2)*(2*continuity+2)] = float[(2*continuity+2)*(2*continuity+2)](1.0000,0.0000,0.0000,0,0,0,1.0000,0.2000,0.0000,0,0,0,1.0000,0.4000,0.0500,0,0,0,0,0,0,1.0000,-0.4000,0.0500,0,0,0,1.0000,-0.2000,-0.0000,0,0,0,1.0000,-0.0000,-0.0000);

vec3[continuity+1] startPointDerivatives(float freedoms[freedom])
{
  uint id = gl_WorkGroupID.x;
  vec3 d0 = freedoms[0] * inBuf.data[2*id].e;
  vec3 dd0 = freedoms[2] * inBuf.data[2*id].e + inBuf.data[2*id].K * freedoms[0] * freedoms[0] * inBuf.data[2*id].n;
  return vec3[continuity+1](inBuf.data[2*id].p, d0, dd0);
}

vec3[continuity+1] endPointDerivatives(float freedoms[freedom])
{
  uint id = gl_WorkGroupID.x;
  vec3 d1 = freedoms[1] * inBuf.data[2*id+1].e;
  vec3 dd1 = freedoms[3] * inBuf.data[2*id+1].e + inBuf.data[2*id+1].K * freedoms[1] * freedoms[1] * inBuf.data[2*id+1].n;
  return vec3[continuity+1](inBuf.data[2*id+1].p, d1, dd1);
}

vec3[point_num] calculateControlPoints(vec3 start[continuity+1], vec3 end[continuity+1])
{
  float pointdata[(2*continuity+2)*3];
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
  float controlpoints[(2*continuity+2)*3] = matmul(pinverse, pointdata);
  vec3 controlPoints[point_num];
  for(int i=0;i < continuity+1;++i)
    controlPoints[i] = vec3(controlpoints[3*i],controlpoints[3*i+1],controlpoints[3*i+2]);
  for(int i=0;i < continuity+1;++i)
    controlPoints[continuity+1+extra_points+i] = vec3(controlpoints[3*(continuity+1)+3*i],controlpoints[3*(continuity+1)+3*i+1],controlpoints[3*(continuity+1)+3*i+2]);
  return controlPoints;
}

#define ITERATIONS 10

float[freedom] optimizeNewton(float start[freedom], bool opt_points, out bool no_move)
{
	uint id = gl_WorkGroupID.x;
  /*uint threadX = gl_LocalInvocationID.x;
  uint threadY = gl_LocalInvocationID.y;
	uint threadZ = gl_LocalInvocationID.z;
  uint threadXY = threadX * 4 + threadY;*/

  float eps = 1e-3;
  float min_dist = length(inBuf.data[2*id].p.xyz - inBuf.data[2*id+1].p.xyz) / 100.0;
  vec3 points[point_num];
  float freedoms[freedom];

  float d[2*continuity];
  float dd[2*continuity*2*continuity];
  float newton_step[2*continuity];
  status _positions[ITERATIONS+1];
  Result positions[ITERATIONS+1];

  for(int i=0;i < freedom;++i) freedoms[i] = start[i];

	for(int i=0;i < ITERATIONS;++i)
	{
    if(!opt_points){
      float integrals[9];
      for(int x=0;x < 2*continuity;++x) for(int y=0;y < 2*continuity;++y){
        for(int z=0;z < 9;++z){
          float freedoms_local[freedom];
          for(int j=0;j < freedom;++j)
          {
            if( j == x ) freedoms_local[j] = freedoms[j]+(int(z)%3-1)*eps;
            else if( j == y ) freedoms_local[j] = freedoms[j]+(int(z)/3-1)*eps;
            else freedoms_local[j] = freedoms[j];
          }
      	  points = calculateControlPoints(startPointDerivatives(freedoms_local), endPointDerivatives(freedoms_local));
          if(infnorm)
      		  integrals[z] = maxvalue(points);
          else
      		  integrals[z] = integral2(points);

          if(x == 0 && y == 0 && z == 4)
          {
            for(int j=0;j < point_num;++j)
      			   positions[i].points[j] = points[j];
      			positions[i].norm = integrals[4];
            _positions[i].freedoms = freedoms;
          }

        }
        d[x] = (integrals[5]-integrals[3])/(2*eps);
        if(x == y) dd[x*2*continuity+y] = (integrals[5] - 2*integrals[4] + integrals[3])/(eps*eps);
        else       dd[x*2*continuity+y] = ((integrals[8]-integrals[6])-(integrals[2]-integrals[0]))/(4*eps*eps);
      }
      newton_step = linsolve(dd,d);
      for(int x=0;x < 2*continuity;++x)
        freedoms[x] -= newton_step[x];

      if(freedoms[0] < min_dist) freedoms[0] = min_dist;
      if(freedoms[1] < min_dist) freedoms[1] = min_dist;

      float n = 0.0;
      for(int j=0;j < 2*continuity;++j)
      {
        n += d[j]*d[j];
      }
      if(n < eps) break;
    }
  }

  points = calculateControlPoints(startPointDerivatives(freedoms), endPointDerivatives(freedoms));
  for(int j=0;j < point_num;++j)
     positions[ITERATIONS].points[j] = points[j];
  if(infnorm)
  	positions[ITERATIONS].norm = maxvalue(points);
  else
    positions[ITERATIONS].norm = integral2(points);
  _positions[ITERATIONS].freedoms = freedoms;

  int minimum = 0;
  //dump.data[0] = positions[0].norm;
  for(int i=1;i<=ITERATIONS;++i)
  {
  	if(!isnan(positions[i].norm) && positions[i].norm < positions[minimum].norm) minimum = i;
    //dump.data[i] = positions[i].norm;
  }
  status _current;
  _current = _positions[minimum];
  //current = positions[minimum];
  float n = 0.f;
  for(int i=0;i < freedom;++i)
  {
    n += pow(start[i]-_current.freedoms[i], 2);
  }
  n = sqrt(n);
  no_move = n < eps;
  return _current.freedoms;
}

//shared float freedoms[freedom];

void main()
{
	uint id = gl_WorkGroupID.x;
  /*uint threadX = gl_LocalInvocationID.x;
  uint threadY = gl_LocalInvocationID.y;
	uint threadZ = gl_LocalInvocationID.z;
  uint threadXY = threadX * 4 + threadY;*/
  vec3 points[point_num];
  float _start[freedom];
  float _current[freedom];
  Result current;

	for(int i=0;i < 2*continuity;++i) _start[i] = length(inBuf.data[2*id].p.xyz - inBuf.data[2*id+1].p.xyz) / 5;
  for(int i=2*continuity;i < freedom;++i) _start[i] = 0;
	points = calculateControlPoints(startPointDerivatives(_start), endPointDerivatives(_start));
  vec3 point_a = points[continuity];
  vec3 point_b = points[point_num-continuity-1];
  for(int i=0; i < extra_points; ++i)
  {
    vec3 p = ((extra_points-i)*point_a + (i+1)*point_b) / (extra_points+1);
    _start[2*continuity+3*i+0] = p.x;
    _start[2*continuity+3*i+1] = p.y;
    _start[2*continuity+3*i+2] = p.z;
  }

  bool finished = false;
  int c = 0;
  while(true)
  {
    _start = _current;

    _current = optimizeNewton(_start, false, finished);

    if(finished) break;

    c++;
  }
  if(_current[0] < length(inBuf.data[2*id].p.xyz - inBuf.data[2*id+1].p.xyz) / 100.0) _current[0] = length(inBuf.data[2*id].p.xyz - inBuf.data[2*id+1].p.xyz) / 100.0;
  if(_current[1] < length(inBuf.data[2*id].p.xyz - inBuf.data[2*id+1].p.xyz) / 100.0) _current[1] = length(inBuf.data[2*id].p.xyz - inBuf.data[2*id+1].p.xyz) / 100.0;
  for(int i=0;i < freedom;++i) dump.data[i] = _current[i];
  dump.data[11] = length(inBuf.data[2*id].p.xyz - inBuf.data[2*id+1].p.xyz) / 100.0;
  points = calculateControlPoints(startPointDerivatives(_current), endPointDerivatives(_current));
  //if(id==1) for(int i=0;i < 12;++i) dump.data[i] = positions[min].points[i/3][i%3];

  for(int j=0;j < point_num;++j)
    outBuf.data[id].points[j] = points[j];
	outBuf.data[id].norm = integral2(points);

	//dump.data[0] = 42;

}
