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
