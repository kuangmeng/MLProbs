#include "PartitionFunction.h"
#include "ExpPartitionFunctionParams.h"
#include "LogPartitionFunctionParams.h"
#include "DataStructures/Sequence.h"
#include "Pairwise/AminoAcidMatrices.h"
#include "Pairwise/NucleotideMatrices.h"
#include "Pairwise/ScoringFactory.h"
#include "Configuration.h"

#include <memory>

using namespace quickprobs;

/// <summary>
/// See declaration for all the details.
/// </summary>
PartitionFunction::PartitionFunction(const Configuration& config) 
{
	auto scoring = ScoringFactory<double>::create(
		config.partition.matrix, 
		config.partition.gapExtend, 
		config.partition.gapOpen);

	if (config.optimisation.useDoublePartition) {
		params = std::make_shared<ExpPartitionFunctionParams<double>>(*scoring, config.partition.temperature);

	} else {
		params = std::make_shared<LogPartitionFunctionParams>(*scoring, config.partition.temperature);
	}
}

/// <summary>
/// See declaration for all the details.
/// </summary>
std::unique_ptr<std::vector<float>> PartitionFunction::computePosteriorMatrix(const Sequence& seq1, const Sequence& seq2)
{
	auto forward = computeForward(seq1, seq2, *params);
	auto posterior = computeReverse(seq1, seq2, *params, *forward);

	return posterior;
}

void PartitionFunction::computePosteriorMatrix(
	const Sequence& seq1, 
	const Sequence& seq2,
	double* forward,
	float* posterior) 
{
	computeForward(seq1, seq2, *params, forward);
	computeReverse(seq1, seq2, *params, forward, posterior);
}

/// <summary>
/// See declaration for all the details.
/// </summary>
std::unique_ptr<std::vector<double>> PartitionFunction::computeForward(
	const Sequence& seq1, 
	const Sequence& seq2, 
	const IPartitionFunctionParams& params)
{
	int seq1Length = seq1.GetLength();
	int seq2Length = seq2.GetLength();
	
	auto result = std::unique_ptr<std::vector<double>>(new std::vector<double>((seq1Length + 1) * (seq2Length + 1), 0));

	computeForward(seq1, seq2, params, result->data());

	return result;
}
	
void PartitionFunction::computeForward(
	const Sequence& seq1, 
	const Sequence& seq2, 
	const IPartitionFunctionParams& params,
	double *forward)
{
	const auto& raw = dynamic_cast<const ExpPartitionFunctionParams<double>&>(params).raw;
	
	int seq1Length = seq1.GetLength();
	int seq2Length = seq2.GetLength();
	int lda = seq2Length + 1;

	const char *s1 = seq1.getData() + 1; // + 1 because of '@' character at the beginning of seq.
	const char *s2 = seq2.getData() + 1;

	auto Zm = forward;
	double* buffer = new double[4 * (seq2Length + 1)];
	std::fill_n(buffer, 4 * (seq2Length + 1), 0);
	double* Ze = buffer + 0 * (seq2Length + 1); // two row vectors
	double* Zf = buffer + 2 * (seq2Length + 1); // two row vectors
	double zz;

	//init z
	Zm[0] = 1.00;
	Zf[0] = Ze[0] = 0;
	Zf[1 * lda + 0] = Zm[0] * raw.termGapOpen;
	Ze[0 * lda + 1] = Zm[0] * raw.termGapOpen;

	//>=2ND COL INIT
	for (int j = 2; j <= seq2Length; j++)
	{
		Ze[0 * lda + j] = Ze[0 * lda + j - 1] * raw.termGapExtend;
	}
	
	//1ST ROW/COL INIT
	for (int i = 1; i <= seq1Length; i++)
	{
		for (int j = 1; j <= seq2Length; j++)
		{
			int Si = s1[i - 1] - 'A';// raw.substIndex[s1[i - 1] - 'A'];
			int Tj = s2[j - 1] - 'A'; //raw.substIndex[s2[j - 1] - 'A'];

			double score = raw.subMatrix[Si * 26 + Tj];
			double open0, extend0, open1, extend1;

			open0 = open1 = raw.gapOpen;
			extend0 = extend1 = raw.gapExt;

			//check to see if one of the 2 sequences or both reach the end
			if (i == seq1Length)
			{
				open0 = raw.termGapOpen;
				extend0 = raw.termGapExtend;
			}

			if (j == seq2Length)
			{
				open1 = raw.termGapOpen;
				extend1 = raw.termGapExtend;
			}

			//z computation using open and extend temp vars
			//open0 is gap open in seq0 and open1 is gap open in seq1
			//entend0 is gap extend in seq0 and extend1 is gap extend in seq1
			Ze[1 * lda + j] = Zm[i * lda + j - 1] * open0 + Ze[1 * lda + j - 1] * extend0;
			Zf[1 * lda + j] = Zm[(i - 1) * lda + j] * open1 + Zf[0 * lda + j] * extend1;
			Zm[i * lda + j] = (Zm[(i - 1) * lda + j - 1] + Ze[0 * lda + j - 1] + Zf[0 * lda + j - 1]) * score;
			
			zz = Zm[i * lda + j] + Ze[1 * lda + j] + Zf[1 * lda + j]; 
		} //end for

		for (int t = 0; t <= seq2Length; t++)
		{
			Ze[0 * lda + t] = Ze[1 * lda + t];
			Ze[1 * lda + t] = 0;

			Zf[0 * lda + t] = Zf[1 * lda + t];
			Zf[1 * lda + t] = 0;
		}

		Zf[1 * lda + 0] = 1;
	}				//end for

	//store the sum of zm zf ze (m,n)s in zm's 0,0 th position
	Zm[0] = zz;
	delete [] buffer;
}

/// <summary>
/// See declaration for all the details.
/// </summary>
std::unique_ptr<std::vector<float>> PartitionFunction::computeReverse(
	const Sequence& seq1, 
	const Sequence& seq2, 
	const IPartitionFunctionParams& params, 
	const std::vector<double> &forward)
{
	int seq1Length = seq1.GetLength();
	int seq2Length = seq2.GetLength();

	auto result = std::unique_ptr<std::vector<float>>(
		new std::vector<float>((seq1Length + 1) * (seq2Length + 1), 0)); // seq1 - vertical, seq2 - horizontal
	std::vector<float>& posterior = *result;

	computeReverse(seq1, seq2, params, forward.data(), result->data());

	return result;
}

void PartitionFunction::computeReverse(
	const Sequence& seq1, 
	const Sequence& seq2, 
	const IPartitionFunctionParams& params, 
	const double* forward,
	float* posterior)
{
	const auto& raw = dynamic_cast<const ExpPartitionFunctionParams<double>&>(params).raw;

	int seq1Length = seq1.GetLength();
	int seq2Length = seq2.GetLength();

	const char *s1 = seq1.getData() + 1;
	const char *s2 = seq2.getData() + 1;

	const double* Zfm = forward;

	double* buffer = new double[6 * (seq2Length + 1)];
	std::fill_n(buffer, 6 * (seq2Length + 1), 0);
	double* Zm = buffer + 0 * (seq2Length + 1); // two row vectors
	double* Ze = buffer + 2 * (seq2Length + 1); // two row vectors
	double* Zf = buffer + 4 * (seq2Length + 1); // two row vectors

	int lda = seq2Length + 1;

	//in Zm
	//let:
	//  Zm(0) be the current row being filled/computed
	//  Zm(1) be the previous row
	Zm[1 * lda + seq2Length] = 1;
	Ze[0 * lda + seq2Length] = Zf[0 * lda + seq2Length] = 0;
	Zf[1 * lda + seq2Length] = Zm[1 * lda + seq2Length] * raw.termGapOpen;
	Ze[0 * lda + seq2Length - 1] = Zm[1 * lda + seq2Length] * raw.termGapOpen;

	//>=2ND COL INIT
	for (int j = seq2Length - 2; j >= 0; j--)
	{
		Ze[0 * lda + j] = Ze[0 * lda + j + 1] * raw.termGapExtend;
	}
	
	for (int i = seq1Length - 1; i >= 0; i--)
	{
		for (int j = seq2Length - 1; j >= 0; j--)
		{
			int Si = s1[i] - 'A';//raw.substIndex[s1[i] - 'A'];
			int Tj = s2[j] - 'A';//raw.substIndex[s2[j] - 'A'];
			double scorez = raw.subMatrix[Si * 26 + Tj];

			//endgaps modification aug 10
			double open0, extend0, open1, extend1;

			open0 = open1 = raw.gapOpen;
			extend0 = extend1 = raw.gapExt;

			//check to see if one of the 2 sequences or both reach the end
			if (i == 0)
			{
				open0 = raw.termGapOpen;
				extend0 = raw.termGapExtend;
			}

			if (j == 0)
			{
				open1 = raw.termGapOpen;
				extend1 = raw.termGapExtend;
			}
		
			//2 ROW zE zF ALGORITHM GOES...:
			//Ze[1][j] =Zm[i][j + 1] * exp(beta * open0) + Ze[1][j + 1] *exp(beta * extend0);
			//Zf[1][j] = Zm[i + 1][j] * exp(beta * open1) + Zf[0][j] * exp(beta * extend1);
			//Zm[i][j] = (Zm[i + 1][j + 1] + Zf[0][j + 1] + Ze[0][j + 1]) * exp(beta * scorez);
			//zz = Zm[0][j] + Zf[1][j] + Ze[1][j];

			//lowmem code for merging probability calculating module
			//Here we make use of Zm as a 2 row matrix

			Zf[1 * lda + j] = Zm[1 * lda + j] * open1 + Zf[0 * lda + j] * extend1;
			Ze[1 * lda + j] = Zm[0 * lda + j + 1] * open0 + Ze[1 * lda + j + 1] * extend0;
			Zm[0 * lda + j] = (Zm[1 * lda + j + 1] + Zf[0 * lda + j + 1] + Ze[0 * lda + j + 1]) * scorez;

			double tempvar = Zfm[(i + 1) * lda + j + 1] * Zm[0 * lda + j];
			//divide P(i,j) i.e. pairwise probability by denominator
			tempvar /= (scorez * Zfm[0]);
			float probability = (float) tempvar;

			//store only noticeable probabilities
			if (probability <= 1 && probability >= 0.001)
			{
				//algorithm goes...
				//validprob[i + 1][j + 1] = probability;
				posterior[(i + 1) * lda + j + 1] = probability;
			}
		}	//end of for

		for (int t = 0; t <= seq2Length; t++)
		{
			Ze[0 * lda + t] = Ze[1 * lda + t];
			Ze[1 * lda + t] = 0;

			Zf[0 * lda + t] = Zf[1 * lda + t];
			Zf[1 * lda + t] = 0;

			Zm[1 * lda + t] = Zm[0 * lda + t];
			Zm[0 * lda + t] = 0;

		}
		Zf[0 * lda + seq2Length] = 1;
	}				//end of for

	delete [] buffer;
	posterior[0] = 0;
}
