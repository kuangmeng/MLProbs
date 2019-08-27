////////////////////////////////////////////////////////////////
// MultiSequence.h
//
// Utilities for reading/writing multiple sequence data.
/////////////////////////////////////////////////////////////////
#pragma once
#include <vector>
#include <memory>
#include <string>
#include <set>

#include "Sequence.h"
#include "ISequenceSet.h"

namespace quickprobs
{
// forward declarations


/////////////////////////////////////////////////////////////////
// MultiSequence
//
// Class for multiple sequence alignment input/output.
/////////////////////////////////////////////////////////////////
class MultiSequence : public ISequenceSet
{
public:
	std::vector<quickprobs::Sequence*> sequences;

	/// Retrieves a sequence from the MultiSequence object.
	Sequence* GetSequence (int i){ return sequences[i]; }

	/// Retrieves a sequence from the MultiSequence object (const version).
	const Sequence* GetSequence (int i) const { return sequences[i]; }

	/// Returns the number of sequences.
	::size_t count () const { return sequences.size(); }

	// Gets the length of the first sequence
	inline ::size_t length() const { return sequences[0]->GetLength(); }

	const ::size_t calculateHash() const;

	/// Default constructor.
	MultiSequence () {}

	MultiSequence (MultiSequence&& other) 
	{
		this->sequences = std::move(other.sequences);
	}

	/// Constructor. Load MFA from a filename.
	MultiSequence (const std::string &filename) { LoadMFA (filename); }

	/////////////////////////////////////////////////////////////////
	// MultiSequence::~MultiSequence()
	//
	// Destructor.  Gets rid of sequence objects contained in the
	// multiple alignment.
	/////////////////////////////////////////////////////////////////
	virtual ~MultiSequence();

	/////////////////////////////////////////////////////////////////
	// MultiSequence::LoadMFA()
	//
	// Load MFA from a filename.
	/////////////////////////////////////////////////////////////////
	void LoadMFA (const std::string &filename, bool stripGaps = false);

	/////////////////////////////////////////////////////////////////
	// MultiSequence::LoadMFA()
	//
	// Load MSF from a FileBuffer object.
	/////////////////////////////////////////////////////////////////
//	void ParseMSF (FileBuffer &infile, string header, bool stripGaps = false);

	/////////////////////////////////////////////////////////////////
	// MultiSequence::AddSequence()
	//
	// Add another sequence to an existing sequence list
	/////////////////////////////////////////////////////////////////
	void AddSequence (Sequence *sequence);

	/////////////////////////////////////////////////////////////////
	// MultiSequence::RemoveSequence()
	//
	// Remove a sequence from the MultiSequence
	/////////////////////////////////////////////////////////////////
	void RemoveSequence (int index);
  
	/////////////////////////////////////////////////////////////////
	// MultiSequence::WriteMFA()
	//
	// Write MFA to the outfile.  Allows the user to specify the
	// number of columns for the output.  Also, useIndices determines
	// whether or not the actual sequence comments will be printed
	// out or whether the artificially assigned sequence labels will
	// be used instead.
	/////////////////////////////////////////////////////////////////
	void WriteMFA (std::ostream &outfile) const;

	/////////////////////////////////////////////////////////////////
	// MultiSequence::GetAnnotationChar()
	//
	// Return CLUSTALW annotation for column.
	/////////////////////////////////////////////////////////////////
	char GetAnnotationChar (std::vector<char> &column);

	/////////////////////////////////////////////////////////////////
	// MultiSequence::WriteALN()
	//
	// Write ALN to the outfile.  Allows the user to specify the
	// number of columns for the output.  
	/////////////////////////////////////////////////////////////////
	void WriteALN (std::ostream &outfile, int numColumns = 60);

	/////////////////////////////////////////////////////////////////
	// MultiSequence::SortByHeader()
	//
	// Organizes the sequences according to their sequence headers
	// in ascending order.
	/////////////////////////////////////////////////////////////////
	void SortByHeader ();

	/////////////////////////////////////////////////////////////////
	// MultiSequence::SortByLabel()
	//
	// Organizes the sequences according to their sequence labels
	// in ascending order.
	/////////////////////////////////////////////////////////////////
	void SortByLabel ();

	/////////////////////////////////////////////////////////////////
	// MultiSequence::SaveOrdering()
	//
	// Relabels sequences so as to preserve the current ordering.
	/////////////////////////////////////////////////////////////////
	void SaveOrdering();

	/////////////////////////////////////////////////////////////////
	// MultiSequence::Project()
	//
	// Given a set of indices, extract all sequences from the current
	// MultiSequence object whose index is included in the set.
	// Then, project the multiple alignments down to the desired
	// subset, and return the projection as a new MultiSequence
	// object.
	/////////////////////////////////////////////////////////////////
	std::unique_ptr<MultiSequence> extractSubset(const std::set<int> &indices) const;
};

};
