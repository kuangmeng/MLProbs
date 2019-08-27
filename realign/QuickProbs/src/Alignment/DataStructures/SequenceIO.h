#pragma once
#include <iostream>

namespace quickprobs
{

// forward declaration 
class MultiSequence;

enum SequenceFormat
{
	FASTA, 
	CLUSTALW,
	PHYLIP
};

class SequenceIO 
{
public:
	/// <summary>
	/// Loads sequence set from specified file in specified format.
	/// </summary>
	/// <param name="filename">Name of input file.</param>
	/// <param name="format">Sequence format.</param>
	/// <param name="set">Output sequence set.</param>
	static void load(std::string filename, SequenceFormat format, MultiSequence& set);

	/// <summary>
	/// Saves sequence set to the specified file in specified format.
	/// </summary>
	/// <param name="filename">Name of input file.</param>
	/// <param name="format">Sequence format.</param>
	/// <param name="set">Output sequence set.</param>
	static void save(std::string filename, SequenceFormat format, MultiSequence& set);

	/// <summary>
	/// Loads set of sequences from Fasta file.
	/// </summary>
	/// <param name="input">Input file.</param>
	/// <param name="sequences">Output set of sequences.</param>
	static void loadFasta(std::istream& input, MultiSequence& sequences);

	/// <summary>
	/// Loads set of sequences from ClustalW file.
	/// </summary>
	/// <param name="input">Input file.</param>
	/// <param name="sequences">Output set of sequences.</param>
	static void loadClustal(std::istream& input, MultiSequence& sequences);
	
	/// <summary>
	/// Loads set of sequences from Phylip file.
	/// </summary>
	/// <param name="input">Input file.</param>
	/// <param name="sequences">Output set of sequences.</param>
	static void loadPhylip(std::istream& input, MultiSequence& sequences);

	static void saveFasta(std::ostream& output, const MultiSequence& alignment);

	static void saveClustal(std::ostream& output, const MultiSequence& alignment);
	
	static void savePhylip(std::ostream& output, const MultiSequence& alignment);

protected:
	static const int MAX_LINE_LENGTH = 10000;

	static const int FASTA_LINE_LENGTH = 60;

	static bool checkAndCorrect(MultiSequence& set);
};

};