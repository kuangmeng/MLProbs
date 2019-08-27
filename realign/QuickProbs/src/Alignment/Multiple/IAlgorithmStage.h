#pragma once
#include "Common/Timer.h"
#include "Common/StatisticsProvider.h"
#include "Configuration.h"

namespace quickprobs 
{

class IAlgorithmStage : public StatisticsProvider
{
public: 
	IAlgorithmStage(std::shared_ptr<Configuration> config) : config(config)
	{
	}

protected:

	std::shared_ptr<Configuration> config;
};

};