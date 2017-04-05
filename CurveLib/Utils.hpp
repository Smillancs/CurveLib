#pragma once

#include <string>
#include <sstream>
#include <cassert>
#include "glm/gtx/string_cast.hpp"

class Exception : public std::exception
{
public:
	Exception(const std::string& str):errorString(str){}
	~Exception() throw() {}
	virtual const char* what() const throw(){return errorString.c_str();}

protected:
	std::string errorString;
};

class UnsupportedDerivativeOrder : public Exception{
public:
	UnsupportedDerivativeOrder(int m) : Exception("This order of derivatives is not supported (" + std::to_string(m) + ")"){}
	~UnsupportedDerivativeOrder() throw() {}
};

inline bool assert_double_equal(double a, double b, bool print = true)
{
	double epsilon = 1e-5;
	bool l = true;
	l = fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
	if( fabs(a)<epsilon && fabs(b)<epsilon ) l = true;
	if(!l && print)
	{
		std::stringstream err;
		err << "Assertion failed: " << a << " != " << b << std::endl;
		throw(Exception(err.str()));
	}
	return l;
}

inline bool assert_dvec3_equal(glm::dvec3 x, glm::dvec3 y, bool print = true)
{
	bool l = true;
	if(!assert_double_equal(x.x, y.x, false)) l = false;
	else if(!assert_double_equal(x.y, y.y, false)) l = false;
	else if(!assert_double_equal(x.z, y.z, false)) l = false;
	if(!l && print)
	{
		std::stringstream err;
		err << "Assertion failed: " << glm::to_string(x) << " != " << glm::to_string(y) << std::endl;
		throw(Exception(err.str()));
	}
	return l;
}

inline bool assert_(bool statement)
{
	if(!statement) throw(Exception("Assertion failed."));
	return statement;
}
