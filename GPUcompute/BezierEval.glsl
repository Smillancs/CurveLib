// Bezier curve evaluation functions in GLSL
// To be included into shaders

int factorial(int n)
{
	int temp = 1;
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
	return binomial(n,i)*pow(t,i)*pow(1.0-t,n-i);
}

vec3 BernsteinEval(float t)
{
	vec3 temp = vec3(0,0,0);
	for(int i = 0; i < pointNum; ++i)
	{
		temp += Bernstein(pointNum-1,i,t)*In[i].pos;
	}
	return temp;
}

vec3 ForwardDiff(int order, int idx)
{
	vec3 delta = vec3(0,0,0);
	int alternate = bool(order%2)?-1:1;
	for (int j = 0; j <= order; ++j)
	{
		delta += binomial(order, j) * float(alternate) * In[idx+j].pos;
		alternate *= -1;
	}
	return delta;
}

vec3 Diff(int order, float t)
{
	vec3 eval = vec3(0,0,0);
	int n = pointNum - order - 1;

	for (int j = 0; j <= n; ++j)
	{
		eval += ForwardDiff(order,j) * Bernstein(n,j,t);
	}
	eval *= factorial(pointNum-1)/float(factorial(n));
	return eval;
}
