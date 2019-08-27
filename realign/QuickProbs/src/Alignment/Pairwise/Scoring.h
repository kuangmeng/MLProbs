#pragma once

#include <exception>
#include <memory>
#include <string>

template <class T>
class Scoring 
{
public:
	T ge;
	T gi;
	const T *substitution;
	std::string alphabet; 
	int symbolsCount;

	T get(int i, int j) const { return substitution[i * symbolsCount + j]; }

	Scoring() : substitution(nullptr), gi(0), ge(0) {}
	
	Scoring(const T *substitution, T ge, T gi, std::string alphabet) 
		: substitution(substitution), gi(gi), ge(ge), alphabet(alphabet), symbolsCount(alphabet.length()) {}

};


/// <summary>
/// Represents system for protein scoring.
/// </summary>
template <class T>
class ProteinScoring : public Scoring<T>
{
public:
	/// <summary>
	/// Initializes substitution matrix with given values.
	/// Gap symbol is also included in matrix.
	/// </summary>
	/// <param name="scores">Scoring matrix.</param>
	ProteinScoring(const T *substitution, T ge, T gi) : Scoring<T>(substitution, ge, gi, "ARNDCQEGHILKMFPSTWYVBZX*") {}
};

/// <summary>
/// Represents system for DNA scoring.
/// </summary>
template <class T>
class DnaScoring : public Scoring<T>
{
public:
	DnaScoring(const T *substitution, T ge, T gi) : Scoring<T>(substitution, ge, gi, "ACGT*") {}
};