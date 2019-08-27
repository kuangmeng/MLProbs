#pragma once
#include "Alignment/DataStructures/SparseMatrixType.h"
#include "Alignment/DataStructures/MultiSequence.h"

#include "ProbabilisticModel.h"
#include "IAlgorithmStage.h"
#include "ConstructionStage.h"
#include "IRefinementObserver.h"

#include <random>
#include <list>

namespace quickprobs 
{

	class RefinementBase : public IAlgorithmStage
	{
	public:
		RefinementBase(
			std::shared_ptr<Configuration> config,
			std::shared_ptr<ConstructionStage> constructor) : IAlgorithmStage(config), constructor(constructor) {}

		virtual std::unique_ptr<MultiSequence> operator()(
			const GuideTree& tree,
			const float *seqsWeights,
			const Array<float>& distances,
			const Array<SparseMatrixType*> &sparseMatrices,
			const ProbabilisticModel &model,
			std::unique_ptr<MultiSequence> alignment);

		void registerObserver(IRefinementObserver* o) { observers.push_back(o); }

	protected:
		std::list<IRefinementObserver*> observers;

		std::shared_ptr<ConstructionStage> constructor;

		int currentIter;

		void notifyObservers(const MultiSequence& alignment, int iteration);

		virtual bool initialise(
			const float *seqsWeights,
			const Array<float>& distances,
			const Array<SparseMatrixType*> &sparseMatrices,
			const ProbabilisticModel &model, 
			const MultiSequence &alignment) { return true; }

		virtual bool finalise() { return true; }

		virtual std::unique_ptr<MultiSequence> refine(
			const GuideTree& tree,
			const float *seqsWeights,
			const Array<float>& distances,
			const Array<SparseMatrixType*> &sparseMatrices,
			const ProbabilisticModel &model, 
			std::unique_ptr<MultiSequence> alignment,
			int depth);

		virtual void split(
			const GuideTree& tree,
			const MultiSequence& alignment,
			std::set<int>& groupOne,
			std::set<int>& groupTwo) = 0;

		virtual bool checkAcceptance(
			const MultiSequence& reference,
			const MultiSequence& candidate);

	};

};