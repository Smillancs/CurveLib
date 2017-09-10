#version 430

layout (local_size_x = 9, local_size_y = 1) in;

const int continuity = 1;
const int point_num = (2*continuity+1)+1;
const int max_length = point_num;

struct ReconstructionData1
{
  vec3 p;
  vec3 e;
};

struct Freedom1
{
  float x1;
};

struct PointData
{
  ReconstructionData1 fix;
  Freedom1 fre;
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

#incl ../Assets/BezierEval.glsl

#incl ../Assets/GeomInvariant.glsl

subroutine float GeomInvariantEval(vec3 data[max_length], int point_num, float t);

subroutine ( GeomInvariantEval ) float curvature(vec3 data[max_length], int point_num, float t) {
    return K(data, point_num, t);
}
subroutine ( GeomInvariantEval ) float curvatureD(vec3 data[max_length], int point_num, float t) {
    return dK(data, point_num, t);
}

subroutine uniform GeomInvariantEval Eval;

// Gauss–Legendre quadrature formula
const float legendre_coeffs[8] = float[8](0.050614, 0.111191, 0.156854, 0.181342, 0.181341, 0.156852, 0.111190, 0.050614);
const float legendre_roots[8]  = float[8](0.019855, 0.101667, 0.237235, 0.408284, 0.591719, 0.762768, 0.898334, 0.980145);

float integral2(vec3 points[max_length])
{
	float s = 0.f;

	for(int i = 0; i < 8; ++i)
	{
		float val = Eval(points, point_num, legendre_roots[i]);
		s += val * val * legendre_coeffs[i];
	}

	return sqrt(s);
}

float[4*3] matmul(float pinverse[4*4], float pointdata[4*3])
{
  float res[4*3];
  for(int i=0;i<4;++i)
    for(int j=0;j<3;++j)
    {
      res[3*i+j] = 0;
      for(int k=0;k<4;++k)
        res[3*i+j] += pinverse[4*i+k]*pointdata[3*k+j];
    }
  return res;
}

const float pinverse[4*4] = float[4*4](1,0,0,0, 1,1/3.0,0,0, 0,0,1,-1/3.0, 0,0,1,0);

vec3[point_num] calculateControlPoints(vec3 start[continuity+1], vec3 end[continuity+1])
{
  float pointdata[4*3];
  for(int i=0;i<2;++i)
  {
    pointdata[3*i+0] = start[i].x;
    pointdata[3*i+1] = start[i].y;
    pointdata[3*i+2] = start[i].z;
  }
  for(int i=0;i<2;++i)
  {
    pointdata[6+3*i+0] = end[i].x;
    pointdata[6+3*i+1] = end[i].y;
    pointdata[6+3*i+2] = end[i].z;
  }
  float controlpoints[4*3] = matmul(pinverse, pointdata);
  vec3 controlPoints[point_num];
  for(int i=0;i<4;++i)
    controlPoints[i] = vec3(controlpoints[3*i],controlpoints[3*i+1],controlpoints[3*i+2]);
  //for(int i=0;i<12;++i) dump.data[i] = controlPoints[i/3][i%3];
  return controlPoints;
}

#define ITERATIONS 10

shared vec3 points[point_num];
shared float integrals[9];
shared float x1_0, x1_1;
shared Result positions[ITERATIONS+1];

void main()
{
	uint id = gl_WorkGroupID.x;
	uint thread = gl_LocalInvocationID.x;
	if(thread == 0)
	{
		x1_0 = 1.f;
		x1_1 = 1.f;
    for(int i=0;i<3;++i) dump.data[3+i] = inBuf.data[2*id].e[i];
	}
	float eps = 1e-3;
  barrier();

	for(int i=0;i<ITERATIONS;++i)
	{
    float x1_0_local = x1_0+(int(thread)%3-1)*eps;
    float x1_1_local = x1_1+(int(thread)/3-1)*eps;
    vec3 d0 = x1_0_local * inBuf.data[2*id].e.xyz;
    vec3 d1 = x1_1_local * inBuf.data[2*id+1].e.xyz;
		vec3 points[point_num] = calculateControlPoints(vec3[continuity+1](inBuf.data[2*id].p.xyz, d0), vec3[continuity+1](inBuf.data[2*id+1].p.xyz, d1));
		integrals[thread] = integral2(points);

		barrier();

		if(thread == 0)

		{
      for(int j=0;j<point_num;++j)
			   positions[i].points[j] = points[j];
			positions[i].norm = integrals[4];

			float d0 = (integrals[5]-integrals[3])/(2*eps);
			float d1 = (integrals[7]-integrals[1])/(2*eps);
			float dd00 = (integrals[5] - 2*integrals[4] + integrals[3])/(eps*eps);
			float dd01 = ((integrals[8]-integrals[6])-(integrals[2]-integrals[0]))/(eps*eps);
			float dd11 = (integrals[7] - 2*integrals[4] + integrals[1])/(eps*eps);

      //for(int i=0;i<9;++i) dump.data[i] = integrals[i];
      /*dump.data[0] = d0;
      dump.data[1] = d1;
      dump.data[2] = dd00;
      dump.data[3] = dd01;
      dump.data[4] = dd01;
      dump.data[5] = dd11;*/

			vec2 diff = vec2(d0,d1);
			mat2 diff2 = mat2(dd00, dd01, dd01, dd11);

			vec2 step = inverse(diff2)*diff;

			x1_0 -= step[0];
			x1_1 -= step[1];

		}
    barrier();
	}


	if(thread == 0)
	{
    vec3 d0 = x1_0 * inBuf.data[2*id].e.xyz;
    vec3 d1 = x1_1 * inBuf.data[2*id+1].e.xyz;
		vec3 points[point_num] = calculateControlPoints(vec3[continuity+1](inBuf.data[2*id].p.xyz, d0), vec3[continuity+1](inBuf.data[2*id+1].p.xyz, d1));
    for(int j=0;j<point_num;++j)
		   positions[ITERATIONS].points[j] = points[j];
		positions[ITERATIONS].norm = integral2(points);

		int min = 0;
		for(int i=1;i<=ITERATIONS;++i)
		{
			if(positions[i].norm < positions[min].norm) min = i;
		}

    for(int j=0;j<point_num;++j)
       outBuf.data[id].points[j] = positions[ITERATIONS].points[j];
		outBuf.data[id].norm = positions[ITERATIONS].norm;

		dump.data[0] = 42;
	}
}
