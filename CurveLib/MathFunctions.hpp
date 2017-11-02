#include <cstdint>
#include <array>
//binomial coeff template metaprograms
//************************************
  #define MAX_BINOMIAL_ORDER 10

  template<int N,int K>
  class allBinomialImpl;
  // this class holds all binomial coeffs up to order N
  template <int N>
  class allBinomial
  {
    union Data
    {
      constexpr Data() : impl(allBinomialImpl<N,N>()) {}
      allBinomialImpl<N,N> impl;
      std::array<int, (N+1)*(N+2)/2> arr;
    };
    const Data data;

  public:
    constexpr int operator()(int n, int k)
    {
      return data.arr[n*(n+1)/2 + k];
    }
  };

  template<int n, int k>
  struct Binomial
  {
    static_assert(n >= k && k >= 0);
    const static int value = (Binomial<n-1,k-1>::value + Binomial<n-1,k>::value);
  };

  template<>
  struct Binomial<0,0>
  {
    const static int value = 1;
  };

  template<int n>
  struct Binomial<n,0>
  {
    const static int value = 1;
  };

  template<int n>
  struct Binomial<n,n>
  {
    const static int value = 1;
  };


  template <int N, int K>
  class allBinomialImpl
  {
    const allBinomialImpl<N,K-1> prev;
    const int value = Binomial<N,K>::value;
  };

  template <int N>
  class allBinomialImpl<N,0>
  {
    const allBinomialImpl<N-1,N-1> prev;
    const int value = Binomial<N,0>::value;
  };

  template<>
  class allBinomialImpl<0,0>
  {
    const int value = Binomial<0,0>::value;
  };

  extern allBinomial<MAX_BINOMIAL_ORDER> binomialCoeff;

  // usage: binomialCoeff(5,2)
//************************************
//end of binomial coeff

#include <functional>
#include "Curve.hpp"

double integrate(std::function<double(Curve&,double)> func, Curve::Ptr c, double a = 0, double b = 1);
