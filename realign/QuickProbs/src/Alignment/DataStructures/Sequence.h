/////////////////////////////////////////////////////////////////
// Sequence.h
//
// Class for reading/manipulating single sequence character data.
/////////////////////////////////////////////////////////////////

// Adam Gudys modifications:
// 06.03.2012: Method GetDataPtr2 added.
#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <memory>

namespace quickprobs
{

/////////////////////////////////////////////////////////////////
// Sequence
//
// Class for storing sequence information.
/////////////////////////////////////////////////////////////////
class Sequence 
{
protected:
	bool isValid;                // a boolean indicating whether the sequence data is valid or not
	std::string header;               // string containing the comment line of the FASTA file
	std::vector<char> *data;      // pointer to character data
	int length;                  // length of the sequence
	int sequenceLabel;           // integer sequence label, typically to indicate the ordering of sequences
								// in a Multi-FASTA file
	int inputLabel;              // position of sequence in original input

	/////////////////////////////////////////////////////////////////
	// Sequence::Sequence()
	//
	// Default constructor.  Does nothing.
	/////////////////////////////////////////////////////////////////
	Sequence () : isValid(false), header(""), data(NULL), length(0), sequenceLabel(0), inputLabel(0) {}

public:

	/////////////////////////////////////////////////////////////////
	// Sequence::Sequence()
	//
	// Constructor.  Builds a sequence from existing data.  Note
	// that the data must use one-based indexing where data[0] should
	// be set to '@'.
	/////////////////////////////////////////////////////////////////
	Sequence (std::vector<char> *data, std::string header, int length, int sequenceLabel, int inputLabel);

	/////////////////////////////////////////////////////////////////
	// Sequence::Sequence()
	//
	// Destructor.  Release allocated memory.
	/////////////////////////////////////////////////////////////////
	virtual ~Sequence ();

	/// Gets sequence header.
	std::string getHeader () const { return header; }

	/// Gets first word of the header.
	std::string getName () const { return std::string(header.substr(0, header.find(" "))); }

	/// Returns the iterator to data associated with this sequence. 
	std::vector<char>::iterator getIterator() { return data->begin(); }

	/// Returns the constant iterator to data associated with this sequence.
	std::vector<char>::const_iterator getIterator() const { return data->begin(); }

	/// Returns the pointer to the sequence.
	char *getData() { return data->data(); }

	/// Returns the constant pointer to the sequence.
	const char *getData() const { return data->data(); }

	/// Returns the character at position i (data is stored with one-based indexing).
	char GetPosition (int i) const {
		assert (i >= 1 && i <= length);
		return (*data)[i];
	}

	/// Sets the sequence label to i.
	void SetLabel (int i) { sequenceLabel = inputLabel = i; }

	/// Sets the sequence sorting label to i.
	void SetSortLabel (int i) { sequenceLabel = i; }

	/// Retrieves the input label.
	int GetLabel () const { return inputLabel; }

	/// Retrieves the sorting label.
	int GetSortLabel () const { return sequenceLabel; }

	/// Checks to see if the sequence successfully loaded.
	bool Fail () const { return !isValid; }

	/// Returns the length of the sequence.
	int GetLength () const { return length; }

	/////////////////////////////////////////////////////////////////
	// Sequence::Clone()
	//
	// Returns a new deep copy of the sequence.
	/////////////////////////////////////////////////////////////////
	Sequence *Clone () const;

	/////////////////////////////////////////////////////////////////
	// Sequence::GetRange()
	//
	// Returns a new sequence object consisting of a range of
	// characters from the current sequence.
	/////////////////////////////////////////////////////////////////
	Sequence *GetRange (int start, int end) const;

	/////////////////////////////////////////////////////////////////
	// Sequence::AddGaps()
	//
	// Given an std::vector<char> containing the skeleton for an
	// alignment and the identity of the current character, this
	// routine will create a new sequence with all necessary gaps added.
	// For instance,
	//    alignment = "XXXBBYYYBBYYXX"
	//    id = 'X'
	// will perform the transformation
	//    "ATGCAGTCA" --> "ATGCC---GT--CA"
	//                    (XXXBBYYYBBYYXX)
	/////////////////////////////////////////////////////////////////
	Sequence *AddGaps (std::vector<char> *alignment, char id) const;

	/////////////////////////////////////////////////////////////////
	// Sequence::GetString()
	//
	// Returns the sequence as a string with gaps removed.
	/////////////////////////////////////////////////////////////////
	std::string GetString ();

	/////////////////////////////////////////////////////////////////
	// Sequence::GetMapping()
	//
	// Returns a std::vector<int> containing the indices of every
	// character in the sequence.  For instance, if the data is
	// "ATGCC---GT--CA", the method returns {1,2,3,4,5,9,10,13,14}.
	/////////////////////////////////////////////////////////////////
	std::unique_ptr<std::vector<int>> getMapping () const;

	void getMapping(int* ret) const;

	/////////////////////////////////////////////////////////////////
	// Sequence::Highlight()
	//
	// Changes all positions with score >= cutoff to upper case and
	// all positions with score < cutoff to lower case.
	/////////////////////////////////////////////////////////////////
	void Highlight (const std::vector<float> &scores, const float cutoff);
	
};

};

