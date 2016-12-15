#ifndef VEC_H
#define VEC_H
#include <vector>
#include <cmath>
#include "all_incl.h"

namespace cu
{

template<class T, size_t s>
class vec
{
public:
	vec() 
	{
		for(size_t i = 0; i < s; ++i)
			push_back(0);
	}

	vec(T x, T y , T z) : points({ x,y,z })
	{}

	vec(T x, T y) : points({ x,y })
	{}

	vec(T x): points({x,0})
	{}

	~vec() {}

	vec & operator= (const vec &)
	{

	}

	inline T at(size_t idx)	{ return points.at(idx); }

	inline T x() { return at(0); }
	inline T y() { return at(1); }
	inline T z() { return at(2); }

	inline void _x(T rx) { at(0) = rx; }
	inline void _y(T ry) { at(1) = ry; }
	inline void _z(T rz) { at(2) = rz; }



	double length()
	{
		double s = 0.0;
		for (size_t i = 0; i < s; ++i)
			s += at(i)*at(i);
		return std::pow(s,1.0f/s);
	}
	
	operator float(){ return 0.0; }
	operator const float(){ return 0.0; }

	template<class T, typename N>
	vec operator*(N num) const
	{
		return vec<T,this.size()>(at(0)*num,at(1)*num);
	}
/*
	template<class T,typename N>
	vec<T, 3> operator*(N num) const
	{
		return vec<T, 3>(at(0)*num, at(1)*num, at(2)*num);
	}*/
	
	vec<T, 2> operator+(const vec<T, 2> &v) const
	{
		vec<T, 2> local = this;
		for (size_t i = 0; i < v.size(); ++i)
		{
			local[i]+=v.at(i);
		}
		return local;
	}


	vec<T, 3> operator+(const vec<T, 3> &v) const
	{
		vec<T, 3> local = this;
		for (size_t i = 0; i < v.size(); ++i)
		{
			local[i]+=v.at(i);
		}
		return local;
	}

	vec<T, 2> operator*(const vec<T, 2> &v) const
	{
		vec<T, 2> local = this;
		for (size_t i = 0; i < v.size(); ++i)
		{
			local[i]*=v.at(i);
		}
		return local;
	}


	vec<T, 3> operator*(const vec<T, 3> &v) const
	{
		vec<T,3> local = this;
		for (size_t i = 0; i < v.size(); ++i)
		{
			local[i]*=v.at(i);
		}
		return local;
	}


	double distance(vec<T, 3> &other)
	{
		assert(other.size() != size());

		double s = 0.0;
		for (size_t i = 0; i < 3; ++i)
			s += (other.at(i) - at(i))*(other.at(i) - at(i));
		return std::pow(s, 1.0f / 3.0f);
	}

	double distance(vec<T, 2> &other)
	{
		assert(other.size() != size());

		double s = 0.0;
		for (size_t i = 0; i < 2; ++i)
			s += (other.at(i) - at(i))*(other.at(i) - at(i));
		return std::pow(s, 1.0f / 2.0f);
	}


private:
	void push_back(T p) { push_back(p); }
	std::vector<T> points;

	size_t size()
	{
		return points.size();
	}
};
}
#endif