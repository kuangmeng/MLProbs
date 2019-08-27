#include <memory>

#include "DataStructures/Sequence.h"
#include "DataStructures/MultiSequence.h"
#include "Pairwise/PairHmmFactory.h"
#include "ProbabilisticModel.h"
#include "ScoreType.h"
#include "Configuration.h"

using namespace quickprobs;

/// <summary>
/// See declaration for all the details.
/// </summary>
ProbabilisticModel::ProbabilisticModel(const Configuration& config)
{
	std::shared_ptr<PairHmm> hmm;

	hmm = PairHmmFactory::create(config.hmm.name, config.hmm.gapOpen, config.hmm.gapExtend);

	this->numStates = hmm->numStates;
	this->numInsertStates = hmm->numInsertStates;

	params.init(*hmm);
	
	// save logarithm from initial distribution
	std::transform(hmm->initDistrib, hmm->initDistrib + hmm->numStates, initialDistribution, [](const float &x)->float {
		return LOG(x);
	});

	// logarithm of partition matrix
	std::transform(hmm->trans, hmm->trans + hmm->numStates * hmm->numStates, &transProb[0][0], [](const float &x)->float {
		return LOG(x);
	});

	float log_e5 = LOG(1e-5);
	float log_e10 = LOG(1e-10);
	std::fill(&insProb[0][0], &insProb[0][0] + 256 * 5, log_e5);
	std::fill(&matchProb[0][0], &matchProb[0][0] + 256 * 256, log_e10);

	// read initial state distribution and transition parameters
	for (int i = 0; i < (int) hmm->alphabet.length(); i++){

		int ii[] = { tolower(hmm->alphabet[i]), toupper(hmm->alphabet[i]) };
		for (int d = 0; d < hmm->numInsertStates; d++) {
			insProb[ii[0]][d] =  insProb[ii[1]][d] = LOG(hmm->emitSingle[i]);
		}

		for (int j = 0; j <= i; j++) {

			int jj[] = { tolower(hmm->alphabet[j]), toupper(hmm->alphabet[j]) };
			matchProb[ii[0]][jj[0]] = matchProb[ii[0]][jj[1]] = matchProb[ii[1]][jj[0]] = matchProb[ii[1]][jj[1]] =
				matchProb[jj[0]][ii[0]] = matchProb[jj[0]][ii[1]] = matchProb[jj[1]][ii[0]] = matchProb[jj[1]][ii[1]] = LOG(hmm->emitPairs[i * hmm->symbolsCount + j]);
		}
	}
}


std::unique_ptr<std::vector<float>> quickprobs::ProbabilisticModel::computeForwardMatrix(
	const Sequence& seq1, 
	const Sequence& seq2,
	float& total) const
{
	int height = seq1.GetLength() + 1;
	int width = seq2.GetLength() + 1;
	std::vector<float> *forwardPtr = new std::vector<float>(numStates * height * width, LOG_ZERO);
	this->computeForwardMatrix(seq1, seq2, forwardPtr->data(), total);
	return std::unique_ptr<std::vector<float>>(forwardPtr);
}

/// <summary>
/// See declaration for all the details.
/// </summary>
void ProbabilisticModel::computeForwardMatrix(
	const Sequence & seq1, 
	const Sequence & seq2,
	float* forward,
	float& total) const
{
	const int seq1Length = seq1.GetLength();
	const int seq2Length = seq2.GetLength();

	// retrieve the points to the beginning of each sequence
	auto iter1 = seq1.getIterator();
	auto iter2 = seq2.getIterator();

	// initialization condition
	forward[0 + numStates * (1 * (seq2Length + 1) + 1)] =
			initialDistribution[0]
					+ matchProb[(unsigned char) iter1[1]][(unsigned char) iter2[1]];

	for (int k = 0; k < numInsertStates; k++) {
		forward[2 * k + 1 + numStates * (1 * (seq2Length + 1) + 0)] =
				initialDistribution[2 * k + 1]
						+ insProb[(unsigned char) iter1[1]][k];
		forward[2 * k + 2 + numStates * (0 * (seq2Length + 1) + 1)] =
				initialDistribution[2 * k + 2]
						+ insProb[(unsigned char) iter2[1]][k];
	}

	// remember offset for each index combination
	int ij = 0;
	int i1j = -seq2Length - 1;
	int ij1 = -1;
	int i1j1 = -seq2Length - 2;

	ij *= numStates;
	i1j *= numStates;
	ij1 *= numStates;
	i1j1 *= numStates;

	// compute forward scores
	for (int i = 0; i <= seq1Length; i++) {
		unsigned char c1 = (i == 0) ? '~' : (unsigned char) iter1[i];
		for (int j = 0; j <= seq2Length; j++) {
			unsigned char c2 = (j == 0) ? '~' : (unsigned char) iter2[j];

			if (i > 1 || j > 1) {
				if (i > 0 && j > 0) {
					forward[0 + ij] = forward[0 + i1j1] + transProb[0][0];
					for (int k = 1; k < numStates; k++)
						LOG_PLUS_EQUALS(forward[0 + ij],
								forward[k + i1j1] + transProb[k][0]);
					forward[0 + ij] += matchProb[c1][c2];
				}
				if (i > 0) {
					for (int k = 0; k < numInsertStates; k++)
						forward[2 * k + 1 + ij] = insProb[c1][k]
								+ LOG_ADD(
										forward[0 + i1j]
												+ transProb[0][2 * k + 1],
										forward[2 * k + 1 + i1j]
												+ transProb[2 * k + 1][2 * k
														+ 1]);
				}
				if (j > 0) {
					for (int k = 0; k < numInsertStates; k++)
						forward[2 * k + 2 + ij] = insProb[c2][k]
								+ LOG_ADD(
										forward[0 + ij1]
												+ transProb[0][2 * k + 2],
										forward[2 * k + 2 + ij1]
												+ transProb[2 * k + 2][2 * k
														+ 2]);
				}
			}

			ij += numStates;
			i1j += numStates;
			ij1 += numStates;
			i1j1 += numStates;
		}
	}
}

/// <summary>
/// See declaration for all the details.
/// </summary>
std::unique_ptr<std::vector<float>> ProbabilisticModel::computeBackwardMatrix(
	const Sequence & seq1, 
	const Sequence & seq2, 
	float& total) const
{
	const int height = seq1.GetLength() + 1;
	const int width = seq2.GetLength() + 1;
	std::vector<float> *backwardPtr = new std::vector<float>(numStates * height * width, LOG_ZERO);
	this->computeBackwardMatrix(seq1, seq2, backwardPtr->data(), total);
	return std::unique_ptr<std::vector<float>>(backwardPtr);
}

void ProbabilisticModel::computeBackwardMatrix(
	const Sequence& seq1, 
	const Sequence& seq2,
	float *backward,
	float& total) const
{
	const int seq1Length = seq1.GetLength();
	const int seq2Length = seq2.GetLength();
	auto iter1 = seq1.getIterator();
	auto iter2 = seq2.getIterator();

	// initialization condition
	for (int k = 0; k < numStates; k++)
		backward[numStates * ((seq1Length + 1) * (seq2Length + 1) - 1)
				+ k] = initialDistribution[k];

	// remember offset for each index combination
	int ij = (seq1Length + 1) * (seq2Length + 1) - 1;
	int i1j = ij + seq2Length + 1;
	int ij1 = ij + 1;
	int i1j1 = ij + seq2Length + 2;

	ij *= numStates;
	i1j *= numStates;
	ij1 *= numStates;
	i1j1 *= numStates;

	// compute backward scores
	for (int i = seq1Length; i >= 0; i--) {
		unsigned char c1 =
				(i == seq1Length) ? '~' : (unsigned char) iter1[i + 1];
		for (int j = seq2Length; j >= 0; j--) {
			unsigned char c2 =
					(j == seq2Length) ? '~' : (unsigned char) iter2[j + 1];

			if (i < seq1Length && j < seq2Length) {
				const float ProbXY = backward[0 + i1j1] + matchProb[c1][c2];
				for (int k = 0; k < numStates; k++)
					LOG_PLUS_EQUALS(backward[k + ij],
							ProbXY + transProb[k][0]);
			}
			if (i < seq1Length) {
				for (int k = 0; k < numInsertStates; k++) {
					LOG_PLUS_EQUALS(backward[0 + ij],
							backward[2 * k + 1 + i1j] + insProb[c1][k]
									+ transProb[0][2 * k + 1]);
					LOG_PLUS_EQUALS(backward[2 * k + 1 + ij],
							backward[2 * k + 1 + i1j] + insProb[c1][k]
									+ transProb[2 * k + 1][2 * k + 1]);
				}
			}
			if (j < seq2Length) {
				for (int k = 0; k < numInsertStates; k++) {
					LOG_PLUS_EQUALS(backward[0 + ij],
							backward[2 * k + 2 + ij1] + insProb[c2][k]
									+ transProb[0][2 * k + 2]);
					LOG_PLUS_EQUALS(backward[2 * k + 2 + ij],
							backward[2 * k + 2 + ij1] + insProb[c2][k]
									+ transProb[2 * k + 2][2 * k + 2]);
				}
			}

			ij -= numStates;
			i1j -= numStates;
			ij1 -= numStates;
			i1j1 -= numStates;
		}
	}
}

/// <summary>
/// See declaration for all the details.
/// </summary>
float ProbabilisticModel::computeTotalProbability(
	int seq1Length, 
	int seq2Length, 
	const float* forward, 
	const float* backward) const
{

	// compute total probability
	float totalForwardProb = LOG_ZERO;
	float totalBackwardProb = LOG_ZERO;
	for (int k = 0; k < numStates; k++) {
		LOG_PLUS_EQUALS(totalForwardProb,
				forward[k
						+ numStates
								* ((seq1Length + 1) * (seq2Length + 1) - 1)]
						+ backward[k
								+ numStates
										* ((seq1Length + 1)
												* (seq2Length + 1) - 1)]);
	}

	totalBackwardProb = forward[0
			+ numStates * (1 * (seq2Length + 1) + 1)]
			+ backward[0 + numStates * (1 * (seq2Length + 1) + 1)];

	for (int k = 0; k < numInsertStates; k++) {
		LOG_PLUS_EQUALS(totalBackwardProb,
				forward[2 * k + 1
						+ numStates * (1 * (seq2Length + 1) + 0)]
						+ backward[2 * k + 1
								+ numStates
										* (1 * (seq2Length + 1) + 0)]);
		LOG_PLUS_EQUALS(totalBackwardProb,
				forward[2 * k + 2
						+ numStates * (0 * (seq2Length + 1) + 1)]
						+ backward[2 * k + 2
								+ numStates
										* (0 * (seq2Length + 1) + 1)]);
	}

	//    cerr << totalForwardProb << " " << totalBackwardProb << endl;

	return (totalForwardProb + totalBackwardProb) / 2;
}


/// <summary>
/// See declaration for all the details.
/// </summary>
std::unique_ptr<std::vector<float>> ProbabilisticModel::computePosteriorMatrix(
	const Sequence & seq1, 
	const Sequence & seq2,
	const std::vector<float> &forward, 
	const std::vector<float> &backward,
	float totalProb) const 
{
	const int seq1Length = seq1.GetLength();
	const int seq2Length = seq2.GetLength();
	std::vector<float> *posteriorPtr = new std::vector<float>((seq1Length + 1) * (seq2Length + 1));
	this->computePosteriorMatrix(seq1, seq2, forward.data(), backward.data(), totalProb, posteriorPtr->data());
	return std::unique_ptr<std::vector<float>>(posteriorPtr);
}


void ProbabilisticModel::computePosteriorMatrix(
	const Sequence & seq1, 
	const Sequence & seq2,
	const float* forward, 
	const float* backward,
	float totalProb,
	float* posterior) const
{
	const int seq1Length = seq1.GetLength();
	const int seq2Length = seq2.GetLength();

	totalProb = computeTotalProbability(seq1Length, seq2Length, forward, backward);

	// compute posterior matrices
	int ij = 0;
	if (totalProb == 0) {
		totalProb = 1.0f;
	}
	
	float* ptr = posterior;

	for (int i = 0; i <= seq1Length; i++) {
		for (int j = 0; j <= seq2Length; j++) {
			*(ptr++) = EXP(
					std::min(LOG_ONE, forward[ij] + backward[ij] - totalProb));
			ij += numStates;
		}
	}

	posterior[0] = 0;

}


/// <summary>
/// See declaration for all the details.
/// </summary>
std::pair<std::vector<char> *, float> ProbabilisticModel::computeAlignment(
	int seq1Length,
	int seq2Length, 
	float *posterior) const 
{
	float *twoRows = new float[(seq2Length + 1) * 2];
	assert(twoRows);
	float *oldRow = twoRows;
	float *newRow = twoRows + seq2Length + 1;

	char *tracebackMatrix = new char[(seq1Length + 1) * (seq2Length + 1)];
	assert(tracebackMatrix);
	char *tracebackPtr = tracebackMatrix;

	float* posteriorPtr = posterior + seq2Length + 1;

	// initialization
	for (int i = 0; i <= seq2Length; i++) {
		oldRow[i] = 0;
		*(tracebackPtr++) = 'L';
	}

	// fill in matrix
	for (int i = 1; i <= seq1Length; i++) {

		// initialize left column
		newRow[0] = 0;
		posteriorPtr++;
		*(tracebackPtr++) = 'U';

		// fill in rest of row
		for (int j = 1; j <= seq2Length; j++) {
			ChooseBestOfThree(*(posteriorPtr++) + oldRow[j - 1],
					newRow[j - 1], oldRow[j], 'D', 'L', 'U', &newRow[j],
					tracebackPtr++);
		}

		// swap rows
		float *temp = oldRow;
		oldRow = newRow;
		newRow = temp;
	}

	// store best score
	float total = oldRow[seq2Length];
	delete[] twoRows;

	// compute traceback
	std::vector<char> *alignment = new std::vector<char>;
	assert(alignment);
	int r = seq1Length, c = seq2Length;
	while (r != 0 || c != 0) {
		char ch = tracebackMatrix[r * (seq2Length + 1) + c];
		switch (ch) {
		case 'L':
			c--;
			alignment->push_back('Y');
			break;
		case 'U':
			r--;
			alignment->push_back('X');
			break;
		case 'D':
			c--;
			r--;
			alignment->push_back('B');
			break;
		default:
			assert(false);
		}
	}

	delete[] tracebackMatrix;

	reverse(alignment->begin(), alignment->end());

	return make_pair(alignment, total);
}

/// <summary>
/// See declaration for all the details.
/// </summary>
std::pair<std::vector<char> *, float> ProbabilisticModel::computeAlignmentWithGapPenalties(
	const MultiSequence & align1, 
	const MultiSequence & align2, 
	const std::vector<float> &posterior,
	int numSeqs1, 
	int numSeqs2, 
	float gapOpenPenalty,
	float gapContinuePenalty) const
{
	int seq1Length = align1.GetSequence(0)->GetLength();
	int seq2Length = align2.GetSequence(0)->GetLength();
	std::vector<std::vector<char>::const_iterator> dataPtrs1(
			align1.count());
	std::vector<std::vector<char>::const_iterator> dataPtrs2(
			align2.count());

	// grab character data
	for (int i = 0; i < align1.count(); i++)
		dataPtrs1[i] = align1.GetSequence(i)->getIterator();
	for (int i = 0; i < align2.count(); i++)
		dataPtrs2[i] = align2.GetSequence(i)->getIterator();

	// the number of active sequences at any given column is defined to be the
	// number of non-gap characters in that column; the number of gap opens at
	// any given column is defined to be the number of gap characters in that
	// column where the previous character in the respective sequence was not
	// a gap
	std::vector<int> numActive1(seq1Length + 1), numGapOpens1(
			seq1Length + 1);
	std::vector<int> numActive2(seq2Length + 1), numGapOpens2(
			seq2Length + 1);

	// compute number of active sequences and gap opens for each group
	for (int i = 0; i < align1.count(); i++) {
		auto dataPtr = align1.GetSequence(i)->getIterator();
		numActive1[0] = numGapOpens1[0] = 0;
		for (int j = 1; j <= seq1Length; j++) {
			if (dataPtr[j] != '-') {
				numActive1[j]++;
				numGapOpens1[j] += (j != 1 && dataPtr[j - 1] != '-');
			}
		}
	}
	for (int i = 0; i < align2.count(); i++) {
		auto dataPtr = align2.GetSequence(i)->getIterator();
		numActive2[0] = numGapOpens2[0] = 0;
		for (int j = 1; j <= seq2Length; j++) {
			if (dataPtr[j] != '-') {
				numActive2[j]++;
				numGapOpens2[j] += (j != 1 && dataPtr[j - 1] != '-');
			}
		}
	}

	Array<float> openingPenalty1(numSeqs1 + 1);
	std::vector<float> continuingPenalty1(numSeqs1 + 1);
	Array<float> openingPenalty2(numSeqs1 + 1);
	std::vector<float> continuingPenalty2(numSeqs2 + 1);

	// precompute penalties
	for (int i = 0; i <= numSeqs1; i++)
		for (int j = 0; j <= numSeqs2; j++)
			openingPenalty1[i][j] = i
					* (gapOpenPenalty * j
							+ gapContinuePenalty * (numSeqs2 - j));
	for (int i = 0; i <= numSeqs1; i++)
		continuingPenalty1[i] = i * gapContinuePenalty * numSeqs2;
	for (int i = 0; i <= numSeqs2; i++)
		for (int j = 0; j <= numSeqs1; j++)
			openingPenalty2[i][j] = i
					* (gapOpenPenalty * j
							+ gapContinuePenalty * (numSeqs1 - j));
	for (int i = 0; i <= numSeqs2; i++)
		continuingPenalty2[i] = i * gapContinuePenalty * numSeqs1;

	float *twoRows = new float[6 * (seq2Length + 1)];
	assert(twoRows);
	float *oldRowMatch = twoRows;
	float *newRowMatch = twoRows + (seq2Length + 1);
	float *oldRowInsertX = twoRows + 2 * (seq2Length + 1);
	float *newRowInsertX = twoRows + 3 * (seq2Length + 1);
	float *oldRowInsertY = twoRows + 4 * (seq2Length + 1);
	float *newRowInsertY = twoRows + 5 * (seq2Length + 1);

	char *tracebackMatrix =
			new char[3 * (seq1Length + 1) * (seq2Length + 1)];
	assert(tracebackMatrix);
	char *tracebackPtr = tracebackMatrix;

	std::vector<float>::const_iterator posteriorPtr = posterior.begin() + seq2Length + 1;

	// initialization
	for (int i = 0; i <= seq2Length; i++) {
		oldRowMatch[i] = oldRowInsertX[i] = (i == 0) ? 0 : LOG_ZERO;
		oldRowInsertY[i] =
				(i == 0) ?
						0 :
						oldRowInsertY[i - 1]
								+ continuingPenalty2[numActive2[i]];
		*(tracebackPtr) = *(tracebackPtr + 1) = *(tracebackPtr + 2) = 'Y';
		tracebackPtr += 3;
	}

	// fill in matrix
	for (int i = 1; i <= seq1Length; i++) {

		// initialize left column
		newRowMatch[0] = newRowInsertY[0] = LOG_ZERO;
		newRowInsertX[0] = oldRowInsertX[0]
				+ continuingPenalty1[numActive1[i]];
		posteriorPtr++;
		*(tracebackPtr) = *(tracebackPtr + 1) = *(tracebackPtr + 2) = 'X';
		tracebackPtr += 3;

		// fill in rest of row
		for (int j = 1; j <= seq2Length; j++) {

			// going to MATCH state
			ChooseBestOfThree(oldRowMatch[j - 1], oldRowInsertX[j - 1],
					oldRowInsertY[j - 1], 'M', 'X', 'Y', &newRowMatch[j],
					tracebackPtr++);
			newRowMatch[j] += *(posteriorPtr++);

			// going to INSERT X state
			ChooseBestOfThree(
					oldRowMatch[j]
							+ openingPenalty1[numActive1[i]][numGapOpens2[j]],
					oldRowInsertX[j] + continuingPenalty1[numActive1[i]],
					oldRowInsertY[j]
							+ openingPenalty1[numActive1[i]][numGapOpens2[j]],
					'M', 'X', 'Y', &newRowInsertX[j], tracebackPtr++);

			// going to INSERT Y state
			ChooseBestOfThree(
					newRowMatch[j - 1]
							+ openingPenalty2[numActive2[j]][numGapOpens1[i]],
					newRowInsertX[j - 1]
							+ openingPenalty2[numActive2[j]][numGapOpens1[i]],
					newRowInsertY[j - 1]
							+ continuingPenalty2[numActive2[j]], 'M', 'X',
					'Y', &newRowInsertY[j], tracebackPtr++);
		}

		// swap rows
		float *temp;
		temp = oldRowMatch;
		oldRowMatch = newRowMatch;
		newRowMatch = temp;
		temp = oldRowInsertX;
		oldRowInsertX = newRowInsertX;
		newRowInsertX = temp;
		temp = oldRowInsertY;
		oldRowInsertY = newRowInsertY;
		newRowInsertY = temp;
	}

	// store best score
	float total;
	char matrix;
	ChooseBestOfThree(oldRowMatch[seq2Length], oldRowInsertX[seq2Length],
			oldRowInsertY[seq2Length], 'M', 'X', 'Y', &total, &matrix);

	delete[] twoRows;

	// compute traceback
	std::vector<char> *alignment = new std::vector<char>;
	assert(alignment);
	int r = seq1Length, c = seq2Length;
	while (r != 0 || c != 0) {

		int offset = (matrix == 'M') ? 0 : (matrix == 'X') ? 1 : 2;
		char ch = tracebackMatrix[(r * (seq2Length + 1) + c) * 3 + offset];
		switch (matrix) {
		case 'Y':
			c--;
			alignment->push_back('Y');
			break;
		case 'X':
			r--;
			alignment->push_back('X');
			break;
		case 'M':
			c--;
			r--;
			alignment->push_back('B');
			break;
		default:
			assert(false);
		}
		matrix = ch;
	}

	delete[] tracebackMatrix;

	reverse(alignment->begin(), alignment->end());

	return make_pair(alignment, 1.0f);
}

/// <summary>
/// See declaration for all the details.
/// </summary>
std::pair<std::vector<char> *, float> ProbabilisticModel::computeViterbiAlignment(
	const Sequence & seq1,
	const Sequence & seq2) const
{
	const int seq1Length = seq1.GetLength();
	const int seq2Length = seq2.GetLength();

	// retrieve the points to the beginning of each sequence
	auto iter1 = seq1.getIterator();
	auto iter2 = seq2.getIterator();

	// create viterbi matrix
	std::vector<float> *viterbiPtr = new std::vector<float>(
			numStates * (seq1Length + 1) * (seq2Length + 1), LOG_ZERO);
	assert(viterbiPtr);
	std::vector<float> &viterbi = *viterbiPtr;

	// create traceback matrix
	std::vector<int> *tracebackPtr = new std::vector<int>(
			numStates * (seq1Length + 1) * (seq2Length + 1), -1);
	assert(tracebackPtr);
	std::vector<int> &traceback = *tracebackPtr;

	// initialization condition
	for (int k = 0; k < numStates; k++)
		viterbi[k] = initialDistribution[k];

	// remember offset for each index combination
	int ij = 0;
	int i1j = -seq2Length - 1;
	int ij1 = -1;
	int i1j1 = -seq2Length - 2;

	ij *= numStates;
	i1j *= numStates;
	ij1 *= numStates;
	i1j1 *= numStates;

	// compute viterbi scores
	for (int i = 0; i <= seq1Length; i++) {
		unsigned char c1 = (i == 0) ? '~' : (unsigned char) iter1[i];
		for (int j = 0; j <= seq2Length; j++) {
			unsigned char c2 = (j == 0) ? '~' : (unsigned char) iter2[j];

			if (i > 0 && j > 0) {
				for (int k = 0; k < numStates; k++) {
					float newVal = viterbi[k + i1j1] + transProb[k][0]
							+ matchProb[c1][c2];
					if (viterbi[0 + ij] < newVal) {
						viterbi[0 + ij] = newVal;
						traceback[0 + ij] = k;
					}
				}
			}
			if (i > 0) {
				for (int k = 0; k < numInsertStates; k++) {
					float valFromMatch = insProb[c1][k] + viterbi[0 + i1j]
							+ transProb[0][2 * k + 1];
					float valFromIns = insProb[c1][k]
							+ viterbi[2 * k + 1 + i1j]
							+ transProb[2 * k + 1][2 * k + 1];
					if (valFromMatch >= valFromIns) {
						viterbi[2 * k + 1 + ij] = valFromMatch;
						traceback[2 * k + 1 + ij] = 0;
					} else {
						viterbi[2 * k + 1 + ij] = valFromIns;
						traceback[2 * k + 1 + ij] = 2 * k + 1;
					}
				}
			}
			if (j > 0) {
				for (int k = 0; k < numInsertStates; k++) {
					float valFromMatch = insProb[c2][k] + viterbi[0 + ij1]
							+ transProb[0][2 * k + 2];
					float valFromIns = insProb[c2][k]
							+ viterbi[2 * k + 2 + ij1]
							+ transProb[2 * k + 2][2 * k + 2];
					if (valFromMatch >= valFromIns) {
						viterbi[2 * k + 2 + ij] = valFromMatch;
						traceback[2 * k + 2 + ij] = 0;
					} else {
						viterbi[2 * k + 2 + ij] = valFromIns;
						traceback[2 * k + 2 + ij] = 2 * k + 2;
					}
				}
			}

			ij += numStates;
			i1j += numStates;
			ij1 += numStates;
			i1j1 += numStates;
		}
	}

	// figure out best terminating cell
	float bestProb = LOG_ZERO;
	int state = -1;
	for (int k = 0; k < numStates; k++) {
		float thisProb =
				viterbi[k
						+ numStates
								* ((seq1Length + 1) * (seq2Length + 1) - 1)]
						+ initialDistribution[k];
		if (bestProb < thisProb) {
			bestProb = thisProb;
			state = k;
		}
	}
	assert(state != -1);

	delete viterbiPtr;

	// compute traceback
	std::vector<char> *alignment = new std::vector<char>;
	assert(alignment);
	int r = seq1Length, c = seq2Length;
	while (r != 0 || c != 0) {
		int newState = traceback[state
				+ numStates * (r * (seq2Length + 1) + c)];

		if (state == 0) {
			c--;
			r--;
			alignment->push_back('B');
		} else if (state % 2 == 1) {
			r--;
			alignment->push_back('X');
		} else {
			c--;
			alignment->push_back('Y');
		}

		state = newState;
	}

	delete tracebackPtr;

	reverse(alignment->begin(), alignment->end());

	return make_pair(alignment, bestProb);
}

/// <summary>
/// See declaration for all the details.
/// </summary>
std::unique_ptr<std::vector<float>> ProbabilisticModel::buildPosterior(
	const MultiSequence & align1, 
	const MultiSequence & align2,
	const Array<SparseMatrixType*>& sparseMatrices,
	float cutoff) const
{
	const int seq1Length = align1.GetSequence(0)->GetLength();
	const int seq2Length = align2.GetSequence(0)->GetLength();

	auto posteriorPtr = std::unique_ptr<std::vector<float>>(
		new std::vector<float>((seq1Length + 1) * (seq2Length + 1), 0));

	std::vector<float> &posterior = *posteriorPtr;
	std::vector<float>::iterator postPtr = posterior.begin();

	// for each s in align1
	for (int i = 0; i < align1.count(); i++) {
		int first = align1.GetSequence(i)->GetLabel();
		auto mapping1 = align1.GetSequence(i)->getMapping();

		// for each t in align2
		for (int j = 0; j < align2.count(); j++) {
			int second = align2.GetSequence(j)->GetLabel();
			auto mapping2 = align2.GetSequence(j)->getMapping();

			if (first < second) {

				// get the associated sparse matrix
				auto matrix = sparseMatrices[first][second];

				for (int ii = 1; ii <= matrix->getSeq1Length(); ii++) {
					auto row = matrix->getRowPtr(ii);
					int base = (*mapping1)[ii] * (seq2Length + 1);
					int rowSize = matrix->getRowSize(ii);

					// add in all relevant values
					for (int jj = 0; jj < rowSize; jj++)
						posterior[base + (*mapping2)[row[jj].getColumn()]] +=
								row[jj].getValue();

					// subtract cutoff 
					for (int jj = 0; jj < matrix->getSeq2Length(); jj++)
						posterior[base + (*mapping2)[jj]] -= cutoff;
				}

			} else {

				// get the associated sparse matrix
				SparseMatrixType *matrix = sparseMatrices[second][first];

				for (int jj = 1; jj <= matrix->getSeq1Length(); jj++) {
					auto row = matrix->getRowPtr(jj);
					int base = (*mapping2)[jj];
					int rowSize = matrix->getRowSize(jj);

					// add in all relevant values
					for (int ii = 0; ii < rowSize; ii++)
						posterior[base
								+ (*mapping1)[row[ii].getColumn()]
										* (seq2Length + 1)] +=
								row[ii].getValue();

					// subtract cutoff 
					for (int ii = 0; ii < matrix->getSeq2Length(); ii++)
						posterior[base + (*mapping1)[ii] * (seq2Length + 1)] -=
								cutoff;
				}

			}
		}
	}

	return posteriorPtr;
}

/// <summary>
/// See declaration for all the details.
/// </summary>
std::unique_ptr<std::vector<float>> ProbabilisticModel::buildPosterior(
	const float* seqsWeights, 
	const MultiSequence & align1,
	const MultiSequence & align2,
	const Array<SparseMatrixType*> &sparseMatrices,
	float cutoff) const 
{
	const int seq1Length = align1.GetSequence(0)->GetLength();
	const int seq2Length = align2.GetSequence(0)->GetLength();

	auto posteriorPtr = std::unique_ptr<std::vector<float>>(
		new std::vector<float>((seq1Length + 1) * (seq2Length + 1), 0));
	
	assert(posteriorPtr);
	std::vector<float> &posterior = *posteriorPtr;
	std::vector<float>::iterator postPtr = posterior.begin();

	//compute the total sum of all weights
	float totalWeights = 0;
	for (int i = 0; i < align1.count(); i++) {
		int first = align1.GetSequence(i)->GetLabel();
		float w1 = seqsWeights[first];
		for (int j = 0; j < align2.count(); j++) {
			int second = align2.GetSequence(j)->GetLabel();
			float w2 = seqsWeights[second];

			totalWeights += w1 * w2;
		}
	}
	// for each s in align1
	for (int i = 0; i < align1.count(); i++) {
		int first = align1.GetSequence(i)->GetLabel();
		// FIXME:
		int w1 = seqsWeights[first];
		auto mapping1 = align1.GetSequence(i)->getMapping();
		// for each t in align2
		for (int j = 0; j < align2.count(); j++) {
			int second = align2.GetSequence(j)->GetLabel();
			float w2 = seqsWeights[second];
			auto mapping2 = align2.GetSequence(j)->getMapping();

			float w = (float) (w1 * w2) / totalWeights;
			if (first < second) {

				// get the associated sparse matrix
				SparseMatrixType *matrix = sparseMatrices[first][second];

				for (int ii = 1; ii <= matrix->getSeq1Length(); ii++) {
					auto row = matrix->getRowPtr(ii);
					int base = (*mapping1)[ii] * (seq2Length + 1);
					int rowSize = matrix->getRowSize(ii);

					// add in all relevant values
					for (int jj = 0; jj < rowSize; jj++)
						posterior[base + (*mapping2)[row[jj].getColumn()]] += w
								* row[jj].getValue();

					// subtract cutoff 
					for (int jj = 0; jj < matrix->getSeq2Length(); jj++)
						posterior[base + (*mapping2)[jj]] -= w * cutoff;
				}

			} else {

				// get the associated sparse matrix
				SparseMatrixType *matrix = sparseMatrices[second][first];

				for (int jj = 1; jj <= matrix->getSeq1Length(); jj++) {
					auto row = matrix->getRowPtr(jj);
					int base = (*mapping2)[jj];
					int rowSize = matrix->getRowSize(jj);

					// add in all relevant values
					for (int ii = 0; ii < rowSize; ii++)
						posterior[base
								+ (*mapping1)[row[ii].getColumn()]
										* (seq2Length + 1)] += w
								* row[ii].getValue();

					// subtract cutoff 
					for (int ii = 0; ii < matrix->getSeq2Length(); ii++)
						posterior[base + (*mapping1)[ii] * (seq2Length + 1)] -=
								w * cutoff;
				}

			}
		}
	}

	return posteriorPtr;
}

