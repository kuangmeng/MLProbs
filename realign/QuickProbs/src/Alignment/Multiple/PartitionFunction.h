#pragma once
#include <memory>
#include <vector>

#include "IPartitionFunctionParams.h"

namespace quickprobs
{
class Sequence;
class Configuration;

/// <summary>
/// This class contains all stuff necessary for computing posterior probability matrix
/// on the basis of partition function. This functionality was originally provided by
/// MSAPartProbs.cpp file.
/// </summary>
class PartitionFunction
{
public:
	std::shared_ptr<IPartitionFunctionParams> params;

	/// <summary>
	/// Performs some initializations.
	/// </summary>
	PartitionFunction(const Configuration& config);

	/// <summary>
	/// Computes posterior matrix on the basis of partition function. 
	/// </summary>
	/// <param name=seq1>First sequence.</param>
	/// <param name=seq1>Second sequence.</param>
	/// <returns>Calculated posterior matrix.</returns>
	virtual std::unique_ptr<std::vector<float>> computePosteriorMatrix(
		const Sequence& seq1, 
		const Sequence& seq2);

	virtual void computePosteriorMatrix(
		const Sequence& seq1, 
		const Sequence& seq2,
		double* forward,
		float* posterior);

	/// <summary>
	/// </summary>
	/// <param name=seq1></param>
	/// <param name=seq2></param>
	/// <param name=info></param>
	/// <returns></returns>
	virtual std::unique_ptr<std::vector<double>> computeForward(
		const Sequence& seq1, 
		const Sequence& seq2, 
		const IPartitionFunctionParams& params);

	virtual void computeForward(
		const Sequence& seq1, 
		const Sequence& seq2, 
		const IPartitionFunctionParams& params,
		double *forward);

	/// <summary>
	/// </summary>
	/// <param name=seq1></param>
	/// <param name=seq2></param>
	/// <param name=info></param>
	/// <param name=forward></param>
	/// <returns></returns>
	virtual std::unique_ptr<std::vector<float>> computeReverse(
		const Sequence& seq1, 
		const Sequence& seq2, 
		const IPartitionFunctionParams& params, 
		const std::vector<double> &forward);

	void computeReverse(
		const Sequence& seq1, 
		const Sequence& seq2, 
		const IPartitionFunctionParams& params, 
		const double* forward,
		float* posterior);

protected:
	double beta;
	double temperature;
};

};