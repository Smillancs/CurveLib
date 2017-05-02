// Bezier curve evaluation functions in GLSL
// To be included into shaders

double factorial(int n)
{
	double temp = 1;
	for (int i = 1; i <= n; ++i)
	{
		temp*=i;
	}
	return temp;
}

double pow(double x, int y)
{
	double ret = 1;
	for(int i = 1; i <= y; ++i)
		ret*=x;
	return ret;
}

double binomial(int n, int k)
{
	return factorial(n)/(factorial(k)*factorial(n-k));
}

double Bernstein(int n, int i, double t)
{
	return binomial(n,i)*pow(t,i)*pow(1.0-t,n-i);
}

dvec3 ForwardDiff(vec3 points[max_length], int order, int idx)
{
	dvec3 delta = vec3(0,0,0);
	int alternate = bool(order%2)?-1:1;
	for (int j = 0; j <= order; ++j)
	{
		delta += binomial(order, j) * float(alternate) * points[idx+j];
		alternate *= -1;
	}
	return delta;
}

dvec3 Diff(vec3 points[max_length], int pointNum, int order, double t)
{
	dvec3 eval = vec3(0,0,0);
	int n = pointNum - order - 1;

	for (int j = 0; j <= n; ++j)
	{
		eval += ForwardDiff(points,order,j) * Bernstein(n,j,t);
	}
	eval *= factorial(pointNum-1)/float(factorial(n));
	return eval;
}

//Nurbs book 21
dvec3 AllBernstein(vec3 points[max_length], int n, double u)
{
	double B[max_length];
	B[0] = 1.0;
	double u1 = 1.0 - u;
	for (int j=1; j<n; ++j)
	{
		double saved = 0.0;
		for (int k=0; k<j; ++k)
		{
			double temp = B[k];
			B[k] = saved+u1*temp;
			saved = u*temp;
		}
		B[j] = saved;
	}
	
	dvec3 ret = vec3(0,0,0);
	for (int i=0; i<n; ++i)
	{
		ret+=B[i]*points[i];
	}
	return ret;
}
