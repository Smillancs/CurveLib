// Geometric invariant calculating functions in GLSL
// To be included into shaders
// Requires the `vec3 Diff(vec3 data[], int order, float t, int pointNum)` nth differentiation function before including (up to 4th order)

dvec3 e(vec3 data[max_length], int point_num, double t)
{
	dvec3 normed = (Diff(data, point_num, 1, t) != dvec3(0,0,0))?
		normalize(Diff(data, point_num, 1, t)):
		dvec3(0,0,0);
	return normed;
}

dvec3 b(vec3 data[max_length], int point_num, double t, int limit = 5)
{
	int i = 2;
	while( cross( Diff(data, point_num, 1, t), Diff(data, point_num, i,t) ) == dvec3(0,0,0) && i < limit ) ++i;
	dvec3 crossed = cross(Diff(data, point_num, 1,t), Diff(data, point_num, i,t));
	dvec3 normed = (crossed != dvec3(0,0,0))?
		normalize(crossed):
		dvec3(0,0,0);
	return normed;
}

dvec3 n(vec3 data[max_length], int point_num, double t)
{
	return cross(b(data, point_num, t), e(data, point_num, t));
}

double frenet_coordinate_e(vec3 data[max_length], int point_num, double t, int order)
{
	return dot(Diff(data, point_num, order, t), e(data, point_num, t));
}

double frenet_coordinate_n(vec3 data[max_length], int point_num, double t, int order)
{
	return dot(Diff(data, point_num, order, t), n(data, point_num, t));
}

double frenet_coordinate_b(vec3 data[max_length], int point_num, double t, int order)
{
	return dot(Diff(data, point_num, order, t), b(data, point_num, t));
}

double K(vec3 data[max_length], int point_num, double t)
{
	double x1 = frenet_coordinate_e(data, point_num, t,1);
	double y2 = frenet_coordinate_n(data, point_num, t,2);
	return y2 / (x1 * x1);
}

double dK(vec3 data[max_length], int point_num, double t)
{
	double x1 = frenet_coordinate_e(data, point_num, t,1);
	double x2 = frenet_coordinate_e(data, point_num, t,2);
	double y2 = frenet_coordinate_n(data, point_num, t,2);
	double y3 = frenet_coordinate_n(data, point_num, t,3);
	return (y3 * x1 - 3 * x2 * y2) / (x1 * x1 * x1 * x1);
}

double ddK(vec3 data[max_length], int point_num, double t)
{
	double x1 = frenet_coordinate_e(data, point_num, t,1);
	double x2 = frenet_coordinate_e(data, point_num, t,2);
	double x3 = frenet_coordinate_e(data, point_num, t,3);
	double y2 = frenet_coordinate_n(data, point_num, t,2);
	double y3 = frenet_coordinate_n(data, point_num, t,3);
	double y4 = frenet_coordinate_n(data, point_num, t,4);
	return (y4 * x1 * x1 - 5 * x2 * y3 * x1 + 12 * x2 * x2 * y2 - 4 * y2 * x3 * x1 - 3 * y2 * y2 * y2) / (x1 * x1 * x1 * x1 * x1 * x1);
}

double T(vec3 data[max_length], int point_num, double t)
{
	double x1 = frenet_coordinate_e(data, point_num, t,1);
	double y2 = frenet_coordinate_n(data, point_num, t,2);
	double z3 = frenet_coordinate_b(data, point_num, t,3);
	return (y2 == 0) ? 0.0 : z3 / (x1 * y2);
}

double dT(vec3 data[max_length], int point_num, double t)
{
	double x1 = frenet_coordinate_e(data, point_num, t,1);
	double y2 = frenet_coordinate_n(data, point_num, t,2);
	double y3 = frenet_coordinate_n(data, point_num, t,3);
	double z3 = frenet_coordinate_b(data, point_num, t,3);
	double z4 = frenet_coordinate_b(data, point_num, t,4);
	return (y2 == 0) ? 0.0 : (z4 * y2 - 2 * z3 * y3) / (x1 * x1 * y2 * y2);
}
