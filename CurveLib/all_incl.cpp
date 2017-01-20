#include "all_incl.h"

int64_t factorial(int64_t n)
{
	long long int out = 1;
	for (size_t i = 1; i <= n; ++i)
	{
		out*=i;
	}
	return out;
}

double binomial(int64_t n, int64_t k)
{
	return factorial(n)/static_cast<double>((factorial(k)*factorial(n-k)));
}