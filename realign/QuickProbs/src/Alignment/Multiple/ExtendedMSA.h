#pragma once
#include "BasicMSA.h"

namespace quickprobs {

class ExtendedMSA : public BasicMSA, public IRefinementObserver {
public:
	static void printWelcome();

	ExtendedMSA() : BasicMSA() {}

	/// <summary>
	/// </summary>
	/// <param name=cl></param>
	/// <returns></returns>
	ExtendedMSA(std::shared_ptr<Configuration> config);

	virtual ~ExtendedMSA() {}

	virtual std::unique_ptr<MultiSequence> doAlign(MultiSequence *sequences);

	virtual void iterationDone(const MultiSequence& alignment, int iteration);

protected:
	std::string version;

	virtual void degenerateDistances(Array<float> &distances);

	virtual void buildDistancesHistogram(const Array<float>& distances);
};

}
