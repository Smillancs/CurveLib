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

const int pointNum = 4;

struct {vec3 pos;} In[pointNum];

#incl ../GPUcompute/BezierEval.glsl

#incl ../GPUcompute/GeomInvariant.glsl

vec2 rotate(vec2 v, float a)
{
	return vec2(v.x*cos(a)-v.y*sin(a), v.y*cos(a)+v.x*sin(a));
}

void calculateControlPoints(float d0, float d1)
{
	uint id = gl_WorkGroupID.x;
	vec2 e0 = rotate(normalize(inBuf.data[id].p1-inBuf.data[id].p0), inBuf.data[id].alpha);
	vec2 e1 = rotate(normalize(inBuf.data[id].p1-inBuf.data[id].p0), inBuf.data[id].beta);
	In[0].pos = vec3(inBuf.data[id].p0, 0);
	In[1].pos = vec3(inBuf.data[id].p0 + d0*e0/(pointNum-1), 0);
	In[2].pos = vec3(inBuf.data[id].p1 - d1*e1/(pointNum-1), 0);
	In[3].pos = vec3(inBuf.data[id].p1, 0);
}

float integral2()
{
	float s = 0;
	float step = 0.1;

	for(float i = step / 2; i < 1.; i += step)
	{
		float val = dK(i);
		s += val * val * step;
	}
	return s;
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
		calculateControlPoints(t0+(int(thread)%3-1)*eps, t1+(int(thread)/3-1)*eps);
		integrals[thread] = integral2();

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
		calculateControlPoints(t0, t1);
		float s = integral2();
		positions[ITERATIONS].t0 = t0;
		positions[ITERATIONS].t1 = t1;
		positions[ITERATIONS].norm = integral2();

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
