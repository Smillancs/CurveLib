#pragma once

#include "Curve.hpp"

#include <vector>

class ExampleHandler
{
	static std::vector<Curve::Ptr> examples;
  static Curve::Ptr random;

	static bool ready;

public:
  static Curve::Ptr getP(int i);
	static Curve& get(int i);

	static size_t size();

  static void newRandom(int deg, int dim);

private:
	static void generate();
};
