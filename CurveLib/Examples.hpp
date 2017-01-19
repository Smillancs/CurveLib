#pragma once

#include "Curve.hpp"

#include <vector>

class ExampleHandler
{
	static std::vector<Curve*> examples;
	
	static bool ready;

public:
	static Curve& get(int i);
	
	static size_t size();
	
private:
	static void generate();
};
