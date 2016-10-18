#include "all_incl.h"

long long int factorial(int n)
{
	long long int out = 1;
	for (size_t i = 1; i < n; ++i)
	{
		out*=i;
	}
	return out;
}

double binomial(int n, int k)
{
	return factorial(n)/((long double)factorial(k)*factorial(n-k));
}