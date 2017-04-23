#version 430

layout (local_size_x = 9, local_size_y = 1) in;

struct Input2D3
{
	vec2 p0;
	vec2 p1;
	float alpha;
	float beta;
};

struct Result
{
	float t0;
	float t1;
	float norm;
};

layout(std430, binding = 0) buffer inputBuffer
{

	Input2D3 data[];

} inBuf;

layout(std430, binding = 1) buffer destBuffer
{

	Result data[];

} outBuf;

layout(std430, binding = 2) buffer debugInfo
{
	float data[100];
} dump;

#define max_length 128

const int pointNum = 4;

#incl ../Assets/BezierEval.glsl

#incl ../Assets/GeomInvariant.glsl

vec2 rotate(vec2 v, float a)
{
	return vec2(v.x*cos(a)-v.y*sin(a), v.y*cos(a)+v.x*sin(a));
}

vec3[max_length] calculateControlPoints(float d0, float d1)
{
	vec3 points[max_length];
	uint id = gl_WorkGroupID.x;
	vec2 e0 = rotate(normalize(inBuf.data[id].p1-inBuf.data[id].p0), inBuf.data[id].alpha);
	vec2 e1 = rotate(normalize(inBuf.data[id].p1-inBuf.data[id].p0), inBuf.data[id].beta);
	points[0] = vec3(inBuf.data[id].p0, 0);
	points[1] = vec3(inBuf.data[id].p0 + d0*e0/(pointNum-1), 0);
	points[2] = vec3(inBuf.data[id].p1 - d1*e1/(pointNum-1), 0);
	points[3] = vec3(inBuf.data[id].p1, 0);
	return points;
}

// Gauss–Legendre quadrature formula
const float legendre_coeffs[8] = float[8](0.050614, 0.111191, 0.156854, 0.181342, 0.181341, 0.156852, 0.111190, 0.050614);
const float legendre_roots[8]  = float[8](0.019855, 0.101667, 0.237235, 0.408284, 0.591719, 0.762768, 0.898334, 0.980145);

float integral2(vec3 points[max_length])
{
	float s = 0.f;

	for(int i = 0; i < 8; ++i)
	{
		float val = dK(points, pointNum, legendre_roots[i]);
		s += val * val * legendre_coeffs[i];
	}

	/*float step = 0.01;
	for(float i = step / 2; i < 1.; i += step)
	{
		float val = dK(points, pointNum, i);
		s += val * val * step;
	}*/
	return sqrt(s);
}

#define ITERATIONS 10

shared float integrals[9];
shared float t0, t1;
shared Result positions[ITERATIONS+1];

void main()
{
	uint id = gl_WorkGroupID.x;
	uint thread = gl_LocalInvocationID.x;
	if(thread == 0)
	{
		t0 = 1.f;
		t1 = 1.f;
	}
	float eps = 1e-3;

	for(int i=0;i<ITERATIONS;++i)
	{
		vec3 points[max_length] = calculateControlPoints(t0+(int(thread)%3-1)*eps, t1+(int(thread)/3-1)*eps);
		integrals[thread] = integral2(points);

		barrier();

		if(thread == 0)

		{

			positions[i].t0 = t0;
			positions[i].t1 = t1;
			positions[i].norm = integrals[4];

			float d0 = (integrals[5]-integrals[3])/(2*eps);
			float d1 = (integrals[7]-integrals[1])/(2*eps);
			float dd00 = (integrals[5] - 2*integrals[4] + integrals[3])/(eps*eps);
			float dd01 = ((integrals[8]-integrals[6])-(integrals[2]-integrals[0]))/(eps*eps);
			float dd11 = (integrals[7] - 2*integrals[4] + integrals[1])/(eps*eps);

			vec2 diff = vec2(d0,d1);
			mat2 diff2 = mat2(dd00, dd01, dd01, dd11);

			vec2 step = inverse(diff2)*diff;

			t0 -= step.x;
			t1 -= step.y;

		}
	}


	if(thread == 0)
	{
		vec3 points[max_length] = calculateControlPoints(t0, t1);
		positions[ITERATIONS].t0 = t0;
		positions[ITERATIONS].t1 = t1;
		positions[ITERATIONS].norm = integral2(points);

		int min = 0;
		for(int i=1;i<=ITERATIONS;++i)
		{
			if(positions[i].norm < positions[min].norm) min = i;
		}

		outBuf.data[id].t0 = positions[min].t0;
		outBuf.data[id].t1 = positions[min].t1;
		outBuf.data[id].norm = positions[min].norm;

		dump.data[0] = 42;
	}
}
