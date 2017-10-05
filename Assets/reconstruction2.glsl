 /* This file is generated from reconstruction.php with generate_reconstruction_glsl.sh
 * For permanent modifications, modify the php file */
#version 430

const int continuity = 2;
const int extra_points = 0;
const int point_num = (2*continuity+1)+1+extra_points;
const int max_length = point_num;
const int freedom = 4;

layout (local_size_x = 4, local_size_y = 4, local_size_z = 9) in;

struct ReconstructionData1
{
  vec3 p;
  vec3 e;
};

struct ReconstructionData2
{
  vec3 p;
  vec3 e;
  vec3 n;
  float K;
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

// Gauss–Legendre quadrature formula
const float legendre_coeffs[20] = float[20](0.00880700357, 0.0203007149, 0.03133602417, 0.04163837079, 0.05096505991, 0.05909726598, 0.06584431922, 0.07104805466, 0.07458649324, 0.07637669357, 0.07637669357, 0.07458649324, 0.07104805466, 0.06584431922, 0.05909726598, 0.05096505991, 0.04163837079, 0.03133602417, 0.0203007149, 0.00880700357);
const float legendre_roots[20]  = float[20](0.003435700407, 0.01801403636, 0.04388278587, 0.08044151409, 0.1268340468, 0.1819731596, 0.244566499, 0.3131469556, 0.3861070744, 0.4617367394, 0.5382632606, 0.6138929256, 0.6868530444, 0.755433501, 0.8180268404, 0.8731659532, 0.9195584859, 0.9561172141, 0.9819859636, 0.9965642996);
// Lobatto-quadrature
const float lobatto_roots[10] = float[10](0, 0.040233045899999986, 0.13061306745, 0.2610375251, 0.41736052115, 0.58263947885, 0.7389624749, 0.86938693255, 0.9597669541, 1);
const float lobatto_coeffs[10] = float[10](0.0111111111, 0.0666529954, 0.112444671, 0.1460213418, 0.1637698806, 0.1637698806, 0.1460213418, 0.112444671, 0.0666529954, 0.0111111111);

float integral2(vec3 points[max_length])
{
	float s = 0.f;

	for(int i = 0; i < 20; ++i)
	{
		float val = Eval(points, point_num, legendre_roots[i]);
		s += val * val * legendre_coeffs[i];
	}

	return sqrt(s);
}

float[point_num*3] matmul(float pinverse[point_num*point_num], float pointdata[point_num*3])
{
  float res[point_num*3];
  for(int i=0;i < point_num;++i)
    for(int j=0;j < 3;++j)
    {
      res[3*i+j] = 0;
      for(int k=0;k < point_num;++k)
        res[3*i+j] += pinverse[point_num*i+k]*pointdata[3*k+j];
    }
  return res;
}

float[freedom] linsolve(float mat[freedom*freedom], float vec[freedom])
{
  for(int k=0;k < freedom;++k)
  {
    for(int i=k+1;i < freedom;++i)
    {
      float f = mat[i*freedom+k] / mat[k*freedom+k];
      for(int j=k+1;j < freedom;++j)
      {
        mat[i*freedom+j] -= mat[k*freedom+j] * f;
      }
      vec[i] -= vec[k]*f;
      mat[i*freedom+k] = 0;
    }
  }

  for(int i=freedom-1;i>=0;--i)
  {
    vec[i] /= mat[i*freedom+i];
    mat[i*freedom+i] = 1;
    for(int j=i-1;j>=0;--j)
    {
      vec[j] -= mat[j*freedom+i] * vec[i];
      mat[j*freedom+i] = 0;
    }
  }
  return vec;
}
const float pinverse[point_num*point_num] = float[point_num*point_num](1.0000,0.0000,0.0000,0,0,0,1.0000,0.2000,0.0000,0,0,0,1.0000,0.4000,0.0500,0,0,0,0,0,0,1.0000,-0.4000,0.0500,0,0,0,1.0000,-0.2000,-0.0000,0,0,0,1.0000,-0.0000,-0.0000);


vec3[point_num] calculateControlPoints(vec3 start[continuity+1], vec3 end[continuity+1])
{
  float pointdata[point_num*3];
  for(int i=0;i < continuity+1;++i)
  {
    pointdata[3*i+0] = start[i].x;
    pointdata[3*i+1] = start[i].y;
    pointdata[3*i+2] = start[i].z;
  }
  for(int i=point_num-continuity-1;i < point_num;++i)
  {
    pointdata[3*i+0] = end[i-(point_num-continuity-1)].x;
    pointdata[3*i+1] = end[i-(point_num-continuity-1)].y;
    pointdata[3*i+2] = end[i-(point_num-continuity-1)].z;
  }
  float controlpoints[point_num*3] = matmul(pinverse, pointdata);
  vec3 controlPoints[point_num];
  for(int i=0;i < point_num;++i)
    controlPoints[i] = vec3(controlpoints[3*i],controlpoints[3*i+1],controlpoints[3*i+2]);
  return controlPoints;
}

#define ITERATIONS 100

vec3 points[point_num];
shared float freedoms[freedom];
shared float integrals[4*4*9];
shared float d[freedom];
shared float dd[freedom*freedom];
shared float step[freedom];
shared Result positions[ITERATIONS+1];

void main()
{
	uint id = gl_WorkGroupID.x;
  uint threadX = gl_LocalInvocationID.x;
  uint threadY = gl_LocalInvocationID.y;
	uint threadZ = gl_LocalInvocationID.z;
  uint threadXY = threadX * 4 + threadY;
	if(threadXY == 0 && threadZ == 0)
	{
		for(int i=0;i < 2*continuity;++i) freedoms[i] = 1;
    for(int i=2*continuity;i < freedom;++i) freedoms[i] = 0;
    /*for(int i=0;i < 3;++i) dump.data[6*id+i] = inBuf.data[2*id].p[i];
    for(int i=0;i < 3;++i) dump.data[6*id+3+i] = inBuf.data[2*id+1].p[i];*/
	}
	float eps = 1e-3;
  float min_dist = length(abs(inBuf.data[2*id].p.xyz - inBuf.data[2*id+1].p.xyz)) / 100;
  barrier();

	for(int i=0;i < ITERATIONS;++i)
	{
    float freedoms_local[freedom];
    for(int j=0;j < freedom;++j)
    {
      if( j == threadX ) freedoms_local[j] = freedoms[j]+(int(threadZ)%3-1)*eps;
      else if( j == threadY ) freedoms_local[j] = freedoms[j]+(int(threadZ)/3-1)*eps;
      else freedoms_local[j] = freedoms[j];
    }
    vec3 d0 = freedoms_local[0] * inBuf.data[2*id].e.xyz;
    vec3 d1 = freedoms_local[1] * inBuf.data[2*id+1].e.xyz;
    vec3 dd0 = freedoms_local[2] * inBuf.data[2*id].e.xyz + inBuf.data[2*id].K * freedoms_local[0] * freedoms_local[0] * inBuf.data[2*id].n;
    vec3 dd1 = freedoms_local[3] * inBuf.data[2*id+1].e.xyz + inBuf.data[2*id+1].K * freedoms_local[1] * freedoms_local[1] * inBuf.data[2*id+1].n;
    points = calculateControlPoints(vec3[continuity+1](inBuf.data[2*id].p.xyz, d0, dd0), vec3[continuity+1](inBuf.data[2*id+1].p.xyz, d1, dd1));
		integrals[threadXY*9+threadZ] = integral2(points);

		barrier();

    if(threadXY == 0 && threadZ == 4)
    {
      for(int j=0;j < point_num;++j)
			   positions[i].points[j] = points[j];
			positions[i].norm = integrals[4];
    }

		if(threadZ == 0)

		{
      d[threadX] = (integrals[threadX*9+5]-integrals[threadX*9+3])/(2*eps);
      if(threadX == threadY)
        dd[threadXY] = (integrals[threadXY*9+5] - 2*integrals[threadXY*9+4] + integrals[threadXY*9+3])/(eps*eps);
      else
        dd[threadXY] = ((integrals[threadXY*9+8]-integrals[threadXY*9+6])-(integrals[threadXY*9+2]-integrals[threadXY*9+0]))/(4*eps*eps);

      //for(int i=0;i < 9;++i) dump.data[i] = integrals[i];
      /*dump.data[0] = d[0];
      dump.data[1] = d[1];
      dump.data[2] = dd[0];
      dump.data[3] = dd[1];
      dump.data[4] = dd[2];
      dump.data[5] = dd[3];*/
      barrier();
      if(threadXY == 0)
  			step = linsolve(dd,d);
      barrier();

      /*dump.data[6] = step[0];
      dump.data[7] = step[1];*/
      if(threadY == 0)
        freedoms[threadX] -= step[threadX];
      barrier();
      if(threadXY == 0)
      {
        if(freedoms[0] < min_dist) freedoms[0] = min_dist;
        if(freedoms[1] < min_dist) freedoms[1] = min_dist;
      }

		}
    barrier();
	}


  if(threadXY == 0 && threadZ == 0)
	{
    vec3 d0 = freedoms[0] * inBuf.data[2*id].e.xyz;
    vec3 d1 = freedoms[1] * inBuf.data[2*id+1].e.xyz;
    vec3 dd0 = freedoms[2] * inBuf.data[2*id].e.xyz + inBuf.data[2*id].K * freedoms[0] * freedoms[0] * inBuf.data[2*id].n;
    vec3 dd1 = freedoms[3] * inBuf.data[2*id+1].e.xyz + inBuf.data[2*id+1].K * freedoms[1] * freedoms[1] * inBuf.data[2*id+1].n;
    points = calculateControlPoints(vec3[continuity+1](inBuf.data[2*id].p.xyz, d0, dd0), vec3[continuity+1](inBuf.data[2*id+1].p.xyz, d1, dd1));
    for(int j=0;j < point_num;++j)
		   positions[ITERATIONS].points[j] = points[j];
		positions[ITERATIONS].norm = integral2(points);

		int min = 0;
		for(int i=1;i<=ITERATIONS;++i)
		{
			if(positions[i].norm < positions[min].norm) min = i;
      //dump.data[i] = positions[i].norm;
		}

    //if(id==1) for(int i=0;i < 12;++i) dump.data[i] = positions[min].points[i/3][i%3];

    for(int j=0;j < point_num;++j)
       outBuf.data[id].points[j] = positions[min].points[j];
		outBuf.data[id].norm = positions[min].norm;

		//dump.data[0] = 42;
	}
}
