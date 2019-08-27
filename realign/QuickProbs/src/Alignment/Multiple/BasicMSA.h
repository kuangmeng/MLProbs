#pragma once
#include <memory>
#include "Common/StatisticsProvider.h"
#include "Common/Array.h"

#include "PosteriorStage.h"
#include "ConsistencyStage.h"
#include "ConstructionStage.h"
#include "RefinementBase.h"

namespace quickprobs {

class BasicMSA : public virtual StatisticsProvider
{
public:
	std::shared_ptr<PosteriorStage> posteriorStage;

	std::shared_ptr<ConsistencyStage> consistencyStage;

	std::shared_ptr<ConstructionStage> constructionStage;

	std::shared_ptr<RefinementBase> refinementStage;

	std::shared_ptr<ProbabilisticModel> model;
	
	BasicMSA(std::shared_ptr<Configuration> config);
	
	virtual ~BasicMSA() {}

	virtual void operator()(std::string inputFile, std::string outputFile);
	
	virtual std::unique_ptr<MultiSequence> doAlign(MultiSequence *sequences);

	virtual void reset() { this->clearStats();  }
	
protected:
	
	std::shared_ptr<Configuration> config;

	BasicMSA();

	void WriteAnnotation(
		MultiSequence *alignment, 
		const Array<SparseMatrixType*> &sparseMatrices);

	int ComputeScore(
		const std::vector<std::pair<int, int>>& active, 
		const Array<SparseMatrixType*> &sparseMatrices);

	void computeDatasetStatistics(const MultiSequence& sequences, const float* weights);

	virtual void printUsage() {};
};

};