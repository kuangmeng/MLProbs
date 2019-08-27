#pragma once
#include <vector>
#include <string>
#include <cmath>

#include "Common/mathex.h"
#include "Hardware/Kernel.h"
#include "Hardware/OpenCl.h"
#include "KernelAlignment/KernelFactory.h"

namespace quickprobs
{


class RelaxationKernelWrapper 
{
public:
	const int stripeCount;

	const int stripeLength;

	const int sparseWidth;
	
	std::shared_ptr<clex::Kernel> obj;

	RelaxationKernelWrapper(
		std::shared_ptr<clex::OpenCL> cl,
		std::string name,
		const std::vector<int>& files, 
		std::vector<std::string> defines,
		int stripeCount,
		int stripeLength,
		int sparseWidth) : stripeCount(stripeCount), stripeLength(stripeLength), sparseWidth(sparseWidth)
	{

		defines.push_back("STRIPE_COUNT=" + std::to_string(stripeCount));
		defines.push_back("STRIPE_LENGTH=" + std::to_string(stripeLength));
		defines.push_back("STRIPE_LENGTH_LOG2=" + std::to_string((int)std::log2(stripeLength)));

		obj = KernelFactory::instance(cl).create(files, name, defines);
	}
};

};