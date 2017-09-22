 /* This file is generated from reconstruction.php with generate_reconstruction_glsl.sh
 * For permanent modifications, modify that file */
#version 430

layout (local_size_x = 9, local_size_y = 1) in;

const int continuity = 2;
const int extra_points = 0;
const int point_num = (2*continuity+1)+1+extra_points;
const int max_length = point_num;
const int freedom = 2*continuity + 3*extra_points;

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

#define ITERATIONS 10

shared vec3 points[point_num];
shared float integrals[9];
shared float freedoms[freedom];
shared Result positions[ITERATIONS+1];

void main()
{
	uint id = gl_WorkGroupID.x;
	uint thread = gl_LocalInvocationID.x;
	if(thread == 0)
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
    uint ternary = thread;
    for(int i=0;i < freedom;++i)
    {
      freedoms_local[i] = freedoms[i]+(int(ternary)%3-1)*eps;
      ternary /= 3;
    }
          vec3 d0 = freedoms_local[0] * inBuf.data[2*id].e.xyz;
      vec3 d1 = freedoms_local[1] * inBuf.data[2*id+1].e.xyz;
      vec3 dd0 = freedoms_local[2] * inBuf.data[2*id].e.xyz + inBuf.data[2*id].K * freedoms_local[0] * freedoms_local[0] * inBuf.data[2*id].n;
      vec3 dd1 = freedoms_local[3] * inBuf.data[2*id+1].e.xyz + inBuf.data[2*id+1].K * freedoms_local[1] * freedoms_local[1] * inBuf.data[2*id+1].n;
      points = calculateControlPoints(vec3[continuity+1](inBuf.data[2*id].p.xyz, d0, dd0), vec3[continuity+1](inBuf.data[2*id+1].p.xyz, d1, dd1));
    		integrals[thread] = integral2(points);

		barrier();

		if(thread == 0)

		{
      for(int j=0;j < point_num;++j)
			   positions[i].points[j] = points[j];
			positions[i].norm = integrals[4];

      float d[freedom];
      float dd[freedom*freedom];

			d[0] = (integrals[5]-integrals[3])/(2*eps);
			d[1] = (integrals[7]-integrals[1])/(2*eps);
			dd[0] = (integrals[5] - 2*integrals[4] + integrals[3])/(eps*eps);
			dd[1] = ((integrals[8]-integrals[6])-(integrals[2]-integrals[0]))/(eps*eps);
      dd[2] = ((integrals[8]-integrals[6])-(integrals[2]-integrals[0]))/(eps*eps);
			dd[3] = (integrals[7] - 2*integrals[4] + integrals[1])/(eps*eps);

      //for(int i=0;i < 9;++i) dump.data[i] = integrals[i];
      /*dump.data[0] = d[0];
      dump.data[1] = d[1];
      dump.data[2] = dd[0];
      dump.data[3] = dd[1];
      dump.data[4] = dd[2];
      dump.data[5] = dd[3];*/

			/*vec2 diff = vec2(d0,d1);
			mat2 diff2 = mat2(dd00, dd01, dd01, dd11);*/

			float step[freedom] = linsolve(dd,d);

      /*dump.data[6] = step[0];
      dump.data[7] = step[1];*/

      for(int i=0;i < freedom;++i) freedoms[i] -= step[i];

      if(freedoms[0] < min_dist) freedoms[0] = min_dist;
      if(freedoms[1] < min_dist) freedoms[1] = min_dist;

		}
    barrier();
	}


	if(thread == 0)
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
		}

    if(id==1) for(int i=0;i < 12;++i) dump.data[i] = positions[min].points[i/3][i%3];

    for(int j=0;j < point_num;++j)
       outBuf.data[id].points[j] = positions[min].points[j];
		outBuf.data[id].norm = positions[min].norm;

		//dump.data[0] = 42;
	}
}
