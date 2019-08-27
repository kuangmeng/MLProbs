#pragma once
#include <vector>
#include <stdexcept>

/// <summary>
/// This interface should be implemented by all classes
/// generating unconstrained pairwise alignment between two strings.
/// </summary>
class IPairwiseAligner
{

public:
	/// <summary>
	/// Calculates pairwise alignment score.
	/// </summary>
	/// <param name="P">Vertical string.</param>
	/// <param name="nP">First string length.</param>
	/// <param name="Q">Second string.</param>
	/// <param name="nQ">Second string length.</param>
	/// <param name="scoring">Scoring system.</param>
	/// <returns>Alignment score.</returns>
	virtual int operator()(
		const char *seq1, 
		const char *seq2,
		int seq1Length,
		int seq2Length) = 0;

	/// <summary>
	/// Performs full pairwise alignment (with backtrack step).
	/// </summary>
	/// <param name="P">First string.</param>
	/// <param name="nP">First string length.</param>
	/// <param name="Q">Second string.</param>
	/// <param name="nQ">Second string length.</param>
	/// <param name="scoring">Scoring system.</param>
	/// <param name="backtrack">Backtrack matrix.</param>
	/// <returns>Alignment score.</returns>
	virtual int operator()(
		const char *seq1, 
		const char *seq2,
		int seq1Length,
		int seq2Length,
		std::vector<int>& scores,
		std::vector<int>& backtrack) 
	{
		std::runtime_error("IPairwiseAligner::operator(): action not implemented.");
		return -1;
	};
};