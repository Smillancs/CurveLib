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
