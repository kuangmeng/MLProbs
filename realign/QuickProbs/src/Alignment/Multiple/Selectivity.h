#pragma once
#include <functional>

namespace quickprobs 
{
	struct Selectivity {
		std::function<float(float,float)> function;
		std::function<float(float,float,float)> filter;
		float filter_a;
		float filter_b;
	};
}