#pragma once
#include "Common/mathex.h"

#include "Alignment/Multiple/Configuration.h"
#include "Hardware/OpenCl.h"
#include "Hardware/Kernel.h"

#include <memory>
#include <map>
#include <limits>

namespace quickprobs {

class KernelSet
{
public:
	enum Identifier {
		HMM_ALL,
		HMM_FORWARD,
		HMM_BACKWARD,
		HMM_COMBINE,
		PARTITION_FORWARD,
		PARTITION_REVERSE,
		FINALIZATION,
		SPARSE
	};

	size_t size() const { return items.size(); }
	
	const std::vector<std::shared_ptr<clex::Kernel>>& getItems() const { return items; }

	KernelSet() {
		items.resize(8);
		elementsPerSymbol.resize(8);
	}

	void fillAll(const Configuration& config, std::shared_ptr<clex::OpenCL> cl);

	void fillAllLong(const Configuration& config, std::shared_ptr<clex::OpenCL> cl);

	clex::Kernel& operator()(int id) {
		return *items[id];
	}

	size_t localMemory(int id, int width) {
		return elementsPerSymbol[id] * mathex::ceilround(width, 32);
	}

	size_t maxWorkgroupSize() {
		size_t s = 0;
		for (auto& kernel : items) {
			if (kernel) {
				s = (kernel->maxWorkgroupSize > s) ? kernel->maxWorkgroupSize : s;
			}
		}
		return s; //fixme
	}

protected:
	std::vector<std::shared_ptr<clex::Kernel>> items;

	std::vector<size_t> elementsPerSymbol;
	
	void insert(int id, std::shared_ptr<clex::Kernel> kernel, size_t elementsPerSymbol) {
		this->items[id] = kernel;
		this->elementsPerSymbol[id] = elementsPerSymbol;
	}
};

}