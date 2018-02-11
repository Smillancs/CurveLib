#pragma once

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
    //static_assert(n >= k && k >= 0);
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

double infnorm(std::function<double(Curve&,double)> func, Curve::Ptr c, double a = 0, double b = 1);


template<int continuity>
std::array<float, (2*continuity+2)*3> matmul(std::array<float,(2*continuity+2)*(2*continuity+2)> pinverse, std::array<float,(2*continuity+2)*3> pointdata)
{
  std::array<float, (2*continuity+2)*3> res;
  for(int i=0;i < (2*continuity+2);++i)
    for(int j=0;j < 3;++j)
    {
      res[3*i+j] = 0;
      for(int k=0;k < (2*continuity+2);++k)
        res[3*i+j] += pinverse[(2*continuity+2)*i+k]*pointdata[3*k+j];
    }
  return res;
}

template<size_t size>
std::array<float,size> linsolve(std::array<float,size*size> mat, std::array<float,size> vec)
{
  for(int k=0;k < size;++k)
  {
    for(int i=k+1;i < size;++i)
    {
      float f = (mat[k*size+k] == 0 ? 0 : mat[i*size+k] / mat[k*size+k]);
      for(int j=k+1;j < size;++j)
      {
        mat[i*size+j] -= mat[k*size+j] * f;
      }
      vec[i] -= vec[k]*f;
      mat[i*size+k] = 0;
    }
  }

  for(int i=(int)size-1;i>=0;--i)
  {
    if(mat[i*size+i] != 0)
    {
      vec[i] /= mat[i*size+i];
      mat[i*size+i] = 1;
    }
    for(int j=i-1;j>=0;--j)
    {
      vec[j] -= mat[j*size+i] * vec[i];
      mat[j*size+i] = 0;
    }
  }
  return vec;
}
