#include "Pairwise/SmithWaterman.h"
#include "Pairwise/NeedlemanWunsch.h"
#include "Pairwise/AminoAcidMatrices.h"
#include "Pairwise/Backtrack.h"

#include "DataStructures/SparseMatrixType.h"
#include "DataStructures/ISequenceSet.h"
#include "DataStructures/Sequence.h"
#include "DataStructures/MultiSequence.h"
#include "MetaPosteriorStage.h"

int indices[] = {0, 20, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16, -1, 19, 17, 22, 18, 21};

using namespace quickprobs;

quickprobs::MetaPosteriorStage::MetaPosteriorStage(
	std::shared_ptr<clex::OpenCL> cl, 
	std::shared_ptr<Configuration> config) : QuickPosteriorStage(cl, config)
{

}

void quickprobs::MetaPosteriorStage::run(
	MultiSequence& sequences, 
	std::vector<std::vector<float>>& distances,
	std::vector<std::vector<SparseMatrixType*>>& matrices)
{
	QuickPosteriorStage::run(sequences, distances, matrices);

	Backtrack bt;
	Blosum50<int> nwSystem(-8, -8);
	NeedlemanWunsch nw(nwSystem.substitution, nwSystem.symbolsCount, nwSystem.gi, nwSystem.ge);
	
	Pam250<int> swSystem;
	SmithWaterman sw(swSystem.substitution, swSystem.symbolsCount, swSystem.gi, swSystem.ge);

	// add additional information
	for (ptrdiff_t i = 0; i < matrices.size(); ++i) {
		for (ptrdiff_t j = i + 1; j < matrices[i].size(); ++j) {
			auto sparseMatrix = matrices[i][j];
			auto dense = sparseMatrix->getPosterior();

			auto seq1 = sequences.GetSequence(i);
			auto seq2 = sequences.GetSequence(j);
			int seq1Length = seq1->GetLength();
			int seq2Length = seq2->GetLength();

			char *p = new char[seq1Length];
			char *q = new char[seq2Length];

			for (int k = 0; k < seq1Length; ++k) {
				p[k] = indices[seq1->GetPosition(k + 1) - 'A'];
			}

			for (int k = 0; k < seq2Length; ++k) {
				q[k] = indices[seq2->GetPosition(k + 1) - 'A'];
			}

			std::vector<int> scores;
			std::vector<int> backtrack;
			std::vector<int> dotMatrix;

			// calculate pid
			nw(p, q, seq1Length, seq2Length, scores, backtrack);
			int matchCount = bt.calculateMatchCount(backtrack, p, q, seq2Length + 1, seq1Length + 1, seq2Length, seq1Length);
			float pid = (float)(matchCount) / std::min(seq1Length, seq2Length);

			float weight = 0.01 * (1 - pid);
		
			//if (pid < 0.2) {
				std::fill(scores.begin(), scores.end(), 0);
				sw(p, q, seq1Length, seq2Length, scores, backtrack);
				bt(backtrack, seq2Length + 1, seq1Length + 1, sw.start_x, sw.start_y, dotMatrix);

				std::vector<float> normalisedScores(scores.size());
				std::transform(dotMatrix.begin(), dotMatrix.end(), normalisedScores.begin(), [](int x)->float {
					return ((float)x); 
				});

				std::transform(dense->begin(), dense->end(), normalisedScores.begin(), dense->begin(), [weight](float x, float y)->float{
					return std::sqrt(std::pow(x,2) * (1 - weight) + y * weight);
				});
			
				matrices[i][j] = new SparseMatrixType(
					sparseMatrix->getSeq1Length(),
					sparseMatrix->getSeq2Length(), 
					dense->data(), 
					config->posteriorCutoff);

				matrices[j][i] = matrices[i][j]->computeTranspose();
		//	}

			delete [] p;
			delete [] q;
		}
	}


}