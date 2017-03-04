// Geometric invariant calculating functions in GLSL
// To be included into shaders
// Requires the `vec3 Diff(int order, float t)` nth differentiation function before including (up to 4th order)

vec3 e(float t)
{
	vec3 normed = (Diff(1,t)!=vec3(0,0,0))?
		normalize(Diff(1,t)):
		vec3(0,0,0);
	return normed;
}

vec3 b(float t, int limit = 5)
{
	int i = 2;
	while( cross( Diff(1,t), Diff(i,t) ) == vec3(0,0,0) && i < limit ) ++i;
	vec3 crossed = cross(Diff(1,t), Diff(i,t));
	vec3 normed = (crossed != vec3(0,0,0))?
		normalize(crossed):
		vec3(0,0,0);
	return normed;
}

vec3 n(float t)
{
	return cross(b(t), e(t));
}

float frenet_coordinate_e(float t, int order)
{
	return dot(Diff(order,t), e(t));
}

float frenet_coordinate_n(float t, int order)
{
	return dot(Diff(order,t), n(t));
}

float frenet_coordinate_b(float t, int order)
{
	return dot(Diff(order,t), b(t));
}

float K(float t)
{
	float x1 = frenet_coordinate_e(t,1);
	float y2 = frenet_coordinate_n(t,2);
	return y2 / (x1 * x1);
}

float dK(float t)
{
	float x1 = frenet_coordinate_e(t,1);
	float x2 = frenet_coordinate_e(t,2);
	float y2 = frenet_coordinate_n(t,2);
	float y3 = frenet_coordinate_n(t,3);
	return (y3 * x1 - 3 * x2 * y2) / (x1 * x1 * x1 * x1);
}

float ddK(float t)
{
	float x1 = frenet_coordinate_e(t,1);
	float x2 = frenet_coordinate_e(t,2);
	float x3 = frenet_coordinate_e(t,3);
	float y2 = frenet_coordinate_n(t,2);
	float y3 = frenet_coordinate_n(t,3);
	float y4 = frenet_coordinate_n(t,4);
	return (y4 * x1 * x1 - 5 * x2 * y3 * x1 + 12 * x2 * x2 * y2 - 4 * y2 * x3 * x1 - 3 * y2 * y2 * y2) / (x1 * x1 * x1 * x1 * x1 * x1);
}

float T(float t)
{
	float x1 = frenet_coordinate_e(t,1);
	float y2 = frenet_coordinate_n(t,2);
	float z3 = frenet_coordinate_b(t,3);
	return (y2 == 0) ? 0.0 : z3 / (x1 * y2);
}

float dT(float t)
{
	float x1 = frenet_coordinate_e(t,1);
	float y2 = frenet_coordinate_n(t,2);
	float y3 = frenet_coordinate_n(t,3);
	float z3 = frenet_coordinate_b(t,3);
	float z4 = frenet_coordinate_b(t,4);
	return (y2 == 0) ? 0.0 : (z4 * y2 - 2 * z3 * y3) / (x1 * x1 * y2 * y2);
}
