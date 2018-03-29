 /* This file is generated from reconstruction.php with generate_reconstruction_glsl.sh
 * For permanent modifications, modify the php file */
#version 430

const int continuity = 1;
const int extra_points = 1;
const int point_num = (2*continuity+1)+1+extra_points;
const int max_length = point_num;
const int freedom = 5;

layout (local_size_x = 1, local_size_y = 1, local_size_z = 1) in;


struct ReconstructionData1
{
  vec3 p;
  vec3 e;
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
	ReconstructionData1 data[]; // will be paired i.e. startpoint, endpoint, startpoint, ...
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

// Bezier curve evaluation functions in GLSL
// To be included into shaders

float factorial(int n)
{
	float temp = 1;
	for (int i = 1; i <= n; ++i)
	{
		temp*=i;
	}
	return temp;
}

float binomial(int n, int k)
{
	return factorial(n)/(factorial(k)*factorial(n-k));
}

float Bernstein(int n, int i, float t)
{
  if(i == 0 && t == 0 || i == n && t == 1) return 1;
	return binomial(n,i)*pow(t,i)*pow(1.0-t,n-i);
}

vec3 BernsteinEval(vec3 points[max_length], float t, int pointNum)
{
	vec3 temp = vec3(0,0,0);
	for(int i = 0; i < pointNum; ++i)
	{
		temp += Bernstein(pointNum-1,i,t)*points[i];
	}
	return temp;
}

vec3 ForwardDiff(vec3 points[max_length], int order, int idx)
{
	vec3 delta = vec3(0,0,0);
	int alternate = bool(order%2)?-1:1;
	for (int j = 0; j <= order; ++j)
	{
		delta += binomial(order, j) * float(alternate) * points[idx+j];
		alternate *= -1;
	}
	return delta;
}

vec3 Diff(vec3 points[max_length], int order, float t, int pointNum)
{
	vec3 eval = vec3(0,0,0);
	int n = pointNum - order - 1;

	for (int j = 0; j <= n; ++j)
	{
		eval += ForwardDiff(points,order,j) * Bernstein(n,j,t);
	}
	eval *= factorial(pointNum-1)/float(factorial(n));
	return eval;
}

//Nurbs book 21
vec3 AllBernstein(vec3 points[max_length], float u, int n)
{
	float B[max_length];
	B[0] = 1.0;
	float u1 = 1.0 - u;
	for (int j=1; j<n; ++j)
	{
		float saved = 0.0;
		for (int k=0; k<j; ++k)
		{
			float temp = B[k];
			B[k] = saved+u1*temp;
			saved = u*temp;
		}
		B[j] = saved;
	}

	vec3 ret = vec3(0,0,0);
	for (int i=0; i<n; ++i)
	{
		ret+=B[i]*points[i];
	}
	return ret;
}
// Geometric invariant calculating functions in GLSL
// To be included into shaders
// Requires the `vec3 Diff(vec3 data[], int order, float t, int pointNum)` nth differentiation function before including (up to 4th order)

vec3 e(vec3 data[max_length], int point_num, float t)
{
	vec3 normed = (Diff(data, 1, t, point_num) != vec3(0,0,0))?
		normalize(Diff(data, 1, t, point_num)):
		vec3(0,0,0);
	return normed;
}

vec3 b(vec3 data[max_length], int point_num, float t, int limit = 5)
{
	int i = 2;
	while( cross( Diff(data, 1,t, point_num), Diff(data, i,t, point_num) ) == vec3(0,0,0) && i < limit ) ++i;
	vec3 crossed = cross(Diff(data, 1,t, point_num), Diff(data, i,t, point_num));
	vec3 normed = (crossed != vec3(0,0,0))?
		normalize(crossed):
		vec3(0,0,0);
	return normed;
}

vec3 n(vec3 data[max_length], int point_num, float t)
{
	return cross(b(data, point_num, t), e(data, point_num, t));
}

float v(vec3 data[max_length], int point_num, float t)
{
  return length(Diff(data, 1, t, point_num));
}

float frenet_coordinate_e(vec3 data[max_length], int point_num, float t, int order)
{
	return dot(Diff(data, order, t, point_num), e(data, point_num, t));
}

float frenet_coordinate_n(vec3 data[max_length], int point_num, float t, int order)
{
	return dot(Diff(data, order, t, point_num), n(data, point_num, t));
}

float frenet_coordinate_b(vec3 data[max_length], int point_num, float t, int order)
{
	return dot(Diff(data, order, t, point_num), b(data, point_num, t));
}

float K(vec3 data[max_length], int point_num, float t)
{
	float x1 = frenet_coordinate_e(data, point_num, t,1);
	float y2 = frenet_coordinate_n(data, point_num, t,2);
	return y2 / (x1 * x1);
}

float dK(vec3 data[max_length], int point_num, float t)
{
	float x1 = frenet_coordinate_e(data, point_num, t,1);
	float x2 = frenet_coordinate_e(data, point_num, t,2);
	float y2 = frenet_coordinate_n(data, point_num, t,2);
	float y3 = frenet_coordinate_n(data, point_num, t,3);
	return (y3 * x1 - 3 * x2 * y2) / (x1 * x1 * x1 * x1);
}

float ddK(vec3 data[max_length], int point_num, float t)
{
	float x1 = frenet_coordinate_e(data, point_num, t,1);
	float x2 = frenet_coordinate_e(data, point_num, t,2);
	float x3 = frenet_coordinate_e(data, point_num, t,3);
	float y2 = frenet_coordinate_n(data, point_num, t,2);
	float y3 = frenet_coordinate_n(data, point_num, t,3);
	float y4 = frenet_coordinate_n(data, point_num, t,4);
	return (y4 * x1 * x1 - 5 * x2 * y3 * x1 + 12 * x2 * x2 * y2 - 4 * y2 * x3 * x1 - 3 * y2 * y2 * y2) / (x1 * x1 * x1 * x1 * x1 * x1);
}

float T(vec3 data[max_length], int point_num, float t)
{
	float x1 = frenet_coordinate_e(data, point_num, t,1);
	float y2 = frenet_coordinate_n(data, point_num, t,2);
	float z3 = frenet_coordinate_b(data, point_num, t,3);
	return (y2 == 0) ? 0.0 : z3 / (x1 * y2);
}

float dT(vec3 data[max_length], int point_num, float t)
{
	float x1 = frenet_coordinate_e(data, point_num, t,1);
	float y2 = frenet_coordinate_n(data, point_num, t,2);
	float y3 = frenet_coordinate_n(data, point_num, t,3);
	float z3 = frenet_coordinate_b(data, point_num, t,3);
	float z4 = frenet_coordinate_b(data, point_num, t,4);
	return (y2 == 0) ? 0.0 : (z4 * y2 - 2 * z3 * y3) / (x1 * x1 * y2 * y2);
}
//#incl ../Assets/BezierEval.glsl

//#incl ../Assets/GeomInvariant.glsl

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
const float lobatto4_roots[4] = float[4](0, 0.2763932023, 0.7236067977, 1);
const float lobatto4_coeffs[4] = float[4](0.0833333333, 0.41666666665, 0.41666666665, 0.0833333333);
const float lobatto8_roots[8] = float[8](0, 0.06412992575, 0.2041499093, 0.39535039105, 0.60464960895, 0.7958500907, 0.93587007425, 1);
const float lobatto8_coeffs[8] = float[8](0.01785714285, 0.10535211355, 0.1705613462, 0.2062293973, 0.2062293973, 0.1705613462, 0.10535211355, 0.01785714285);
const float lobatto16_roots[16] = float[16](0, 0.0152159769, 0.0503997335, 0.1039958541, 0.1738056486, 0.2569702891, 0.35008476555, 0.44933686325, 0.55066313675, 0.64991523445, 0.7430297109, 0.8261943514, 0.8960041459, 0.9496002665, 0.9847840231, 1);
const float lobatto16_coeffs[16] = float[16](0.00416666665, 0.0290149465, 0.04469684865, 0.06212769105, 0.0770134904, 0.08874595665, 0.0968450119, 0.10097915405, 0.10097915405, 0.0968450119, 0.08874595665, 0.0770134904, 0.06212769105, 0.04469684865, 0.0290149465, 0.00416666665);
const float lobatto32_roots[32] = float[32](0, 0.00369553305, 0.0123526548, 0.0258575808, 0.0440750305, 0.066823762, 0.09387763415, 0.12496775305, 0.1597851222, 0.19798370645, 0.2391838686, 0.2829761414, 0.32892529675, 0.3765746706, 0.4254507016, 0.4750676375, 0.5249323625, 0.5745492984, 0.6234253294, 0.67107470325, 0.7170238586, 0.7608161314, 0.80201629355, 0.8402148778, 0.87503224695, 0.90612236585, 0.933176238, 0.9559249695, 0.9741424192, 0.9876473452, 0.99630446695, 1);
const float lobatto32_coeffs[32] = float[32](0.0010080645, 0.00619905325, 0.0110997764, 0.0158875677, 0.02051710075, 0.02494263565, 0.0291202486, 0.0330084386, 0.0365685698, 0.0397652628, 0.04256674895, 0.04494518645, 0.04687693775, 0.04834280445, 0.04932821825, 0.04982338575, 0.04982338575, 0.04932821825, 0.04834280445, 0.04687693775, 0.04494518645, 0.04256674895, 0.0397652628, 0.0365685698, 0.0330084386, 0.0291202486, 0.02494263565, 0.02051710075, 0.0158875677, 0.0110997764, 0.00619905325, 0.0010080645);
const float lobatto10_roots[10] = float[10](0, 0.0402330459, 0.13061306745, 0.2610375251, 0.41736052115, 0.58263947885, 0.7389624749, 0.86938693255, 0.9597669541, 1);
const float lobatto10_coeffs[10] = float[10](0.0111111111, 0.0666529954, 0.112444671, 0.1460213418, 0.1637698806, 0.1637698806, 0.1460213418, 0.112444671, 0.0666529954, 0.0111111111);

float lobatto_root(uint deg, uint num)
{
  if(deg == 4) return lobatto4_roots[num];
  else if(deg == 8) return lobatto8_roots[num];
  else if(deg == 16) return lobatto16_roots[num];
  else if(deg == 32) return lobatto32_roots[num];
}
float lobatto_coeff(uint deg, uint num)
{
  if(deg == 4) return lobatto4_coeffs[num];
  else if(deg == 8) return lobatto8_coeffs[num];
  else if(deg == 16) return lobatto16_coeffs[num];
  else if(deg == 32) return lobatto32_coeffs[num];
}

float integral2(vec3 points[max_length], uint deg)
{
	float s = 0.f;

	for(int i = 0; i < deg; ++i)
	{
		float val = Eval(points, point_num, lobatto_root(deg, i));
    float vel = v(points, point_num, lobatto_root(deg, i));
		s += val * val * vel * vel * lobatto_coeff(deg, i);
	}

	return sqrt(s);
}
float integral2_uniform(vec3 points[max_length], uint num)
{
	float s = 0.0;

  for(int i = 0; i < num; ++i)
  {
    float root = i/float(num);
    float val = Eval(points, point_num, i/float(num));
    float vel = v(points, point_num, i/float(num));
    s += val * val * vel * vel / num;
  }
	return sqrt(s);
}

float integrate(vec3 points[max_length])
{
  float s = 0.f;
  float sprev = -1.f;
  float eps = 1e-3;
  uint deg = 4;
  while(deg <= 8 && abs(s-sprev) > eps)
  {
    sprev = s;
    s = integral2(points, deg);
    deg *= 2;
  }
  return s;
}


float maxvalue(vec3 points[max_length], uint num)
{
  float mx = 0.f;

  for(int i = 1; i < num; ++i)
  {
    float val = abs(Eval(points, point_num, i/float(num)));
    if(val > mx) mx = val;
  }

  return mx;
}

float infinitenorm(vec3 points[max_length])
{
  float s = 0.f;
  float sprev = -1.f;
  float eps = 1e-3;
  uint num = 100;
  while(num <= 10000 && abs(s-sprev) > eps)
  {
    sprev = s;
    s = maxvalue(points, num);
    num *= 2;
  }
  return s;
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

float[2] linsolve(float mat[4], float vec[2])
{
  for(int k=0;k < 2;++k)
  {
    for(int i=k+1;i < 2;++i)
    {
      float f = mat[i*2+k] / mat[k*2+k];
      for(int j=k+1;j < 2;++j)
      {
        mat[i*2+j] -= mat[k*2+j] * f;
      }
      vec[i] -= vec[k]*f;
      mat[i*2+k] = 0;
    }
  }

  for(int i=2-1;i>=0;--i)
  {
    vec[i] /= mat[i*2+i];
    mat[i*2+i] = 1;
    for(int j=i-1;j>=0;--j)
    {
      vec[j] -= mat[j*2+i] * vec[i];
      mat[j*2+i] = 0;
    }
  }
  return vec;
}
float[3] linsolve(float mat[9], float vec[3])
{
  for(int k=0;k < 3;++k)
  {
    for(int i=k+1;i < 3;++i)
    {
      float f = mat[i*3+k] / mat[k*3+k];
      for(int j=k+1;j < 3;++j)
      {
        mat[i*3+j] -= mat[k*3+j] * f;
      }
      vec[i] -= vec[k]*f;
      mat[i*3+k] = 0;
    }
  }

  for(int i=3-1;i>=0;--i)
  {
    vec[i] /= mat[i*3+i];
    mat[i*3+i] = 1;
    for(int j=i-1;j>=0;--j)
    {
      vec[j] -= mat[j*3+i] * vec[i];
      mat[j*3+i] = 0;
    }
  }
  return vec;
}

const float pinverse[(2*continuity+2)*(2*continuity+2)] = float[(2*continuity+2)*(2*continuity+2)](1,0,0,0, 1,1/3.0,0,0, 0,0,1,-1/3.0, 0,0,1,0);

vec3[continuity+1] startPointDerivatives(float freedoms[freedom])
{
  uint id = gl_WorkGroupID.x;
  vec3 d0 = freedoms[0] * inBuf.data[2*id].e;
	return vec3[continuity+1](inBuf.data[2*id].p, d0);
}

vec3[continuity+1] endPointDerivatives(float freedoms[freedom])
{
  uint id = gl_WorkGroupID.x;
  vec3 d1 = freedoms[1] * inBuf.data[2*id+1].e;
	return vec3[continuity+1](inBuf.data[2*id+1].p, d1);
}

vec3[point_num] calculateControlPoints(vec3 start[continuity+1], vec3 end[continuity+1], vec3 extras[extra_points])
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
  for(int i=0;i < extra_points;++i)
    controlPoints[continuity+1+i] = extras[i];
  for(int i=0;i < continuity+1;++i)
    controlPoints[continuity+1+extra_points+i] = vec3(controlpoints[3*(continuity+1)+3*i],controlpoints[3*(continuity+1)+3*i+1],controlpoints[3*(continuity+1)+3*i+2]);
  return controlPoints;
}

#define ITERATIONS 5

float[freedom] optimizeNewton(float start[freedom], bool opt_points, out bool no_move)
{
	uint id = gl_WorkGroupID.x;
  /*uint threadX = gl_LocalInvocationID.x;
  uint threadY = gl_LocalInvocationID.y;
	uint threadZ = gl_LocalInvocationID.z;
  uint threadXY = threadX * 5 + threadY;*/

  float eps = 1e-3;
  float min_dist = 0.1;//length(inBuf.data[2*id].p.xyz - inBuf.data[2*id+1].p.xyz) / 100.0;
  vec3 points[point_num];
  float freedoms[freedom];

  float d[2*continuity];
  float dd[2*continuity*2*continuity];
  float newton_step[2*continuity];
  float _d[3*extra_points];
  float _dd[3*extra_points*3*extra_points];
  float _newton_step[3*extra_points];
  status _positions[ITERATIONS+1];
  Result positions[ITERATIONS+1];

  for(int i=0;i < freedom;++i) freedoms[i] = start[i];

	for(int i=0;i < ITERATIONS;++i)
	{
    if(!opt_points){
      vec3 extras[extra_points];
      for(int j=0;j < extra_points;++j)
        extras[j] = vec3(freedoms[2*continuity+3*j+0],freedoms[2*continuity+3*j+1],freedoms[2*continuity+3*j+2]);
      float integrals[9];
      for(int x=0;x < 2*continuity;++x) for(int y=0;y <= x;++y){
        for(int z=0;z < 9;++z){
          if(z == 1 || z == 7) continue;
          float freedoms_local[freedom];
          for(int j=0;j < freedom;++j)
          {
            if( j == x ) freedoms_local[j] = freedoms[j]+(int(z)%3-1)*eps;
            else if( j == y ) freedoms_local[j] = freedoms[j]+(int(z)/3-1)*eps;
            else freedoms_local[j] = freedoms[j];
          }
      	  points = calculateControlPoints(startPointDerivatives(freedoms_local), endPointDerivatives(freedoms_local), extras);
          if(infnorm)
      		  integrals[z] = infinitenorm(points);
          else
      		  integrals[z] = integrate(points);

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
        else{
          dd[x*2*continuity+y] = ((integrals[8]-integrals[6])-(integrals[2]-integrals[0]))/(4*eps*eps);
          dd[y*2*continuity+x] = ((integrals[8]-integrals[6])-(integrals[2]-integrals[0]))/(4*eps*eps);
        }
      }
      newton_step = linsolve(dd,d);
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
      float integrals[9];
      for(int x=2*continuity;x < freedom;++x) for(int y=2*continuity;y <= x;++y){
        for(int z=0;z < 9;++z){
          if(z == 1 || z == 7) continue;
          float freedoms_local[freedom];
          for(int j=0;j < freedom;++j)
          {
            if( j == x ) freedoms_local[j] = freedoms[j]+(int(z)%3-1)*eps;
            else if( j == y ) freedoms_local[j] = freedoms[j]+(int(z)/3-1)*eps;
            else freedoms_local[j] = freedoms[j];
          }
          vec3 extras[extra_points];
          for(int j=0;j < extra_points;++j)
            extras[j] = vec3(freedoms_local[2*continuity+3*j+0],freedoms_local[2*continuity+3*j+1],freedoms_local[2*continuity+3*j+2]);
          points = calculateControlPoints(startPointDerivatives(freedoms_local), endPointDerivatives(freedoms_local), extras);
          if(infnorm)
      		  integrals[z] = infinitenorm(points);
          else
      		  integrals[z] = integrate(points);

          if(x == 2*continuity && y == 2*continuity && z == 4)
          {
            for(int j=0;j < point_num;++j)
      			   positions[i].points[j] = points[j];
      			positions[i].norm = integrals[4];
            _positions[i].freedoms = freedoms;
          }
        }
        _d[x] = (integrals[5]-integrals[3])/(2*eps);
        if(x == y) _dd[x*2*continuity+y] = (integrals[5] - 2*integrals[4] + integrals[3])/(eps*eps);
        else{
          _dd[x*2*continuity+y] = ((integrals[8]-integrals[6])-(integrals[2]-integrals[0]))/(4*eps*eps);
          _dd[y*2*continuity+x] = ((integrals[8]-integrals[6])-(integrals[2]-integrals[0]))/(4*eps*eps);
        }
      }
      _newton_step = linsolve(_dd,_d);
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

  vec3 extras[extra_points];
  for(int j=0;j < extra_points;++j)
    extras[j] = vec3(freedoms[2*continuity+3*j+0],freedoms[2*continuity+3*j+1],freedoms[2*continuity+3*j+2]);
  points = calculateControlPoints(startPointDerivatives(freedoms), endPointDerivatives(freedoms), extras);
  for(int j=0;j < point_num;++j)
     positions[ITERATIONS].points[j] = points[j];
  if(infnorm)
  	positions[ITERATIONS].norm = infinitenorm(points);
  else
    positions[ITERATIONS].norm = integrate(points);
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
  uint threadXY = threadX * 5 + threadY;*/
  vec3 points[point_num];
  float _start[freedom];
  float _current[freedom];
  Result current;

	for(int i=0;i < 2*continuity;++i) _start[i] = length(inBuf.data[2*id].p.xyz - inBuf.data[2*id+1].p.xyz) / 5;
  for(int i=2*continuity;i < freedom;++i) _start[i] = 0;
  vec3 extras[extra_points];
  for(int j=0;j < extra_points;++j)
    extras[j] = vec3(_start[2*continuity+3*j+0],_start[2*continuity+3*j+1],_start[2*continuity+3*j+2]);
	points = calculateControlPoints(startPointDerivatives(_start), endPointDerivatives(_start), extras);
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
  _current = optimizeNewton(_start, true, finished);
  while(true)
  {
    _start = _current;

    _current = optimizeNewton(_start, false, finished);

    if(finished) break;

    _start = _current;

    _current = optimizeNewton(_start, true, finished);

    if(finished || c >= 2) break;
    c++;
  }
  if(_current[0] < length(inBuf.data[2*id].p.xyz - inBuf.data[2*id+1].p.xyz) / 100.0) _current[0] = length(inBuf.data[2*id].p.xyz - inBuf.data[2*id+1].p.xyz) / 100.0;
  if(_current[1] < length(inBuf.data[2*id].p.xyz - inBuf.data[2*id+1].p.xyz) / 100.0) _current[1] = length(inBuf.data[2*id].p.xyz - inBuf.data[2*id+1].p.xyz) / 100.0;
  for(int i=0;i < freedom;++i) dump.data[i] = _current[i];
  dump.data[11] = length(inBuf.data[2*id].p.xyz - inBuf.data[2*id+1].p.xyz) / 100.0;
  for(int j=0;j < extra_points;++j)
    extras[j] = vec3(_current[2*continuity+3*j+0],_current[2*continuity+3*j+1],_current[2*continuity+3*j+2]);
  points = calculateControlPoints(startPointDerivatives(_current), endPointDerivatives(_current), extras);
  //if(id==1) for(int i=0;i < 12;++i) dump.data[i] = positions[min].points[i/3][i%3];

  for(int j=0;j < point_num;++j)
    outBuf.data[id].points[j] = points[j];
	outBuf.data[id].norm = integrate(points);

	//dump.data[0] = 42;

}
