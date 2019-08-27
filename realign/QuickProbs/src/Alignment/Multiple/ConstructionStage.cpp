#include "ConstructionStage.h"
#include "ParallelProbabilisticModel.h"
#include "DataStructures/Sequence.h"
#include "DataStructures/MultiSequence.h"

#include "Common/Log.h"

using namespace std;
using namespace quickprobs;

std::unique_ptr<quickprobs::MultiSequence> ConstructionStage::operator()(
	const float *seqsWeights,
	const Array<float>& distances,
	GuideTree *tree, 
	const quickprobs::MultiSequence& sequences, 
	const Array<SparseMatrixType*> &sparseMatrices, 
	const quickprobs::ProbabilisticModel &model)
{
	TIMER_CREATE(timer);
	TIMER_START(timer);
	
	auto out = run(seqsWeights, distances, tree, sequences, sparseMatrices, model);

	TIMER_STOP(timer);

	return out;
}

std::unique_ptr<quickprobs::MultiSequence> ConstructionStage::run(
	const float *seqsWeights,
	const Array<float>& distances,
	GuideTree* tree, 
	const quickprobs::MultiSequence& sequences, 
	const Array<SparseMatrixType*> &sparseMatrices, 
	const quickprobs::ProbabilisticModel &model)
{
	auto alignment = processTree(seqsWeights, distances, *tree->getRoot(), sequences, sparseMatrices, model);

	std::vector<int> oldOrdering;
	if (config->io.enableAlignOrder){
		for (int i = 0; i < alignment->count(); i++) {
			oldOrdering.push_back(alignment->GetSequence(i)->GetSortLabel());
		}
		alignment->SaveOrdering();
		config->io.enableAlignOrder = false;
	}
	
	return alignment;
}


std::unique_ptr<quickprobs::MultiSequence> ConstructionStage::processTree (
	const float *seqsWeights,
	const Array<float>& distances,
	const Node& tree, 
	const quickprobs::MultiSequence& sequences,
	const Array<SparseMatrixType*> &sparseMatrices,
	const quickprobs::ProbabilisticModel &model)
{
	std::unique_ptr<MultiSequence> result;

	// check if this is a node of the alignment tree
	//if (tree->GetSequenceLabel() == -1){
	
	if(tree.type == INTERNAL_NODE) {
		auto alignLeft = processTree(seqsWeights, distances, *(tree.left), sequences, sparseMatrices, model);
		auto alignRight = processTree(seqsWeights, distances, *(tree.right), sequences, sparseMatrices, model);
		assert (alignLeft);
		assert (alignRight);

		result = alignAlignments(seqsWeights, distances, *alignLeft, *alignRight, sparseMatrices, model, nullptr);
		assert (result);

	} else {
		// otherwise, this is a leaf of the alignment tree
		result = std::unique_ptr<MultiSequence>(new MultiSequence()); 
		//result->AddSequence (sequences->GetSequence(tree->GetSequenceLabel())->Clone());
		result->AddSequence (sequences.GetSequence(tree.idx)->Clone());
	}

	return result;
}

std::unique_ptr<quickprobs::MultiSequence> ConstructionStage::alignAlignments (
	const float *seqsWeights,
	const Array<float>& distances,
	const quickprobs::MultiSequence& align1, 
	const quickprobs::MultiSequence& align2,
	const Array<SparseMatrixType*> &sparseMatrices,
	const quickprobs::ProbabilisticModel &model,
	float * buffer)
{
	const ParallelProbabilisticModel &ppm = dynamic_cast<const ParallelProbabilisticModel&>(model);
	
	TIMER_CREATE(timer);
	TIMER_START(timer);	

	std::unique_ptr<std::vector<float>> uptr;

	if (buffer == nullptr) {
		uptr = ppm.buildPosterior(seqsWeights, distances, align1, align2, sparseMatrices, selectivity);
		buffer = uptr->data();
	}  else {
		ppm.buildPosterior(seqsWeights, distances, align1, align2, sparseMatrices, selectivity, buffer);
	}
	TIMER_STOP(timer);
	timeMatrixAddition += timer.seconds();

	//perform alignment
	TIMER_START(timer);
	auto alignment = ppm.computeAlignment(align1.GetSequence(0)->GetLength(), align2.GetSequence(0)->GetLength(), buffer);
	TIMER_STOP(timer);
	timeDynamicProgramming += timer.seconds();

	// now build final alignment
	auto result = std::unique_ptr<MultiSequence>(new MultiSequence());
	for (int i = 0; i < align1.count(); i++)
		result->AddSequence (align1.GetSequence(i)->AddGaps(alignment.first, 'X'));
	for (int i = 0; i < align2.count(); i++)
		result->AddSequence (align2.GetSequence(i)->AddGaps(alignment.first, 'Y'));
	
	if (!config->io.enableAlignOrder)
		result->SortByLabel();

	// free temporary alignment
	delete alignment.first;

	return result;
}