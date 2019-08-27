////////////////////////////////////////////////////////////////
// MultiSequence.h
//
// Utilities for reading/writing multiple sequence data.
/////////////////////////////////////////////////////////////////
#include <cctype>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <set>

#include "Sequence.h"
#include "MultiSequence.h"
#include "SequenceIO.h"

using namespace std;
using namespace quickprobs;


MultiSequence::~MultiSequence()
{
	// free all sequences
	for (std::vector<Sequence *>::iterator iter = sequences.begin(); iter != sequences.end(); ++iter){
		assert (*iter);
		delete *iter;
		*iter = NULL;
	}
}

void MultiSequence::LoadMFA (const std::string &filename, bool stripGaps)
{
	SequenceIO::load(filename, FASTA, *this);
}
/*
void MultiSequence::ParseMSF (FileBuffer &infile, string header, bool stripGaps)
{
	std::vector<std::vector<char> *> seqData;
	std::vector<string> seqNames;
	std::vector<int> seqLengths;

	istringstream in;
	bool valid = true;
	bool missingHeader = false;
	bool clustalW = false;

	// read until data starts
	while (!infile.eof() && header.find ("..", 0) == string::npos){
		if (header.find ("CLUSTAL", 0) == 0 || header.find ("MSAPROBS", 0) == 0){
			clustalW = true;
			break;
		}
		infile.GetLine (header);
		if (header.find ("//", 0) != string::npos){
			missingHeader = true;
			break;
		}
	}

	// read until end-of-file
	while (valid){
		infile.GetLine (header);
		if (infile.eof()) break;

		string word;
		in.clear();
		in.str(header);

		// check if there's anything on this line
		if (in >> word){

			// clustalw name parsing
			if (clustalW){
				if (!isspace(header[0]) && find (seqNames.begin(), seqNames.end(), word) == seqNames.end()){
					seqNames.push_back (word);
					seqData.push_back (new std::vector<char>());
					seqLengths.push_back (0);
					seqData[(int) seqData.size() - 1]->push_back ('@');
				}	  
			}

			// look for new sequence label
			if (word == string ("Name:")){
				if (in >> word){
					seqNames.push_back (word);
					seqData.push_back (new std::vector<char>());
					seqLengths.push_back (0);
					seqData[(int) seqData.size() - 1]->push_back ('@');
				}
				else
					valid = false;
			}

			// check if this is sequence data
			else if (find (seqNames.begin(), seqNames.end(), word) != seqNames.end()){
				int index = find (seqNames.begin(), seqNames.end(), word) - seqNames.begin();

				// read all remaining characters on the line
				char ch;
				while (in >> ch){
					if (isspace (ch)) continue;
					if (ch >= 'a' && ch <= 'z') ch = ch - 'a' + 'A';
					if (ch == '.') ch = '-';
					if (stripGaps && ch == '-') continue;
					if (!((ch >= 'A' && ch <= 'Z') || ch == '*' || ch == '-')){
						cerr << "ERROR: Unknown character encountered: " << ch << endl;
						exit (1);
					}

					// everything's ok so far, so just store this character.
					seqData[index]->push_back (ch);
					seqLengths[index]++;
				}
			}
			else if (missingHeader){
				seqNames.push_back (word);
				seqData.push_back (new std::vector<char>());
				seqLengths.push_back (0);
				seqData[(int) seqData.size() - 1]->push_back ('@');

				int index = (int) seqNames.size() - 1;

				// read all remaining characters on the line
				char ch;
				while (in >> ch){
					if (isspace (ch)) continue;
					if (ch >= 'a' && ch <= 'z') ch = ch - 'a' + 'A';
					if (ch == '.') ch = '-';
					if (stripGaps && ch == '-') continue;
					if (!((ch >= 'A' && ch <= 'Z') || ch == '*' || ch == '-')){
						cerr << "ERROR: Unknown character encountered: " << ch << endl;
						exit (1);
					}

					// everything's ok so far, so just store this character.
					seqData[index]->push_back (ch);
					seqLengths[index]++;
				}
			}
		}
	}

	// check for errors
	if (seqNames.size() == 0){
		cerr << "ERROR: No sequences read!" << endl;
		exit (1);
	}

	assert (!sequences);
	sequences = new std::vector<Sequence *>;
	for (int i = 0; i < (int) seqNames.size(); i++){
		if (seqLengths[i] == 0){
			cerr << "ERROR: Sequence of zero length!" << endl;
			exit (1);
		}
		Sequence *seq = new Sequence (seqData[i], seqNames[i], seqLengths[i], i, i);
		sequences->push_back (seq);
	}
}
*/
void MultiSequence::AddSequence (Sequence *sequence){
	assert (sequence);
	assert (!sequence->Fail());

	// add sequence
	sequences.push_back (sequence);
}

void MultiSequence::RemoveSequence (int index){
	assert (index >= 0 && index < (int) sequences.size());
	delete sequences[index];

	sequences.erase (sequences.begin() + index);
}

void MultiSequence::WriteMFA (std::ostream &outfile) const
{
	SequenceIO::saveFasta(outfile, *this);
}

char MultiSequence::GetAnnotationChar (std::vector<char> &column){
	std::vector<int> counts (256, 0);
	int allChars = (int) column.size();

	for (int i = 0; i < allChars; i++){
		counts[(unsigned char) toupper(column[i])]++;
	}

	allChars -= counts[(unsigned char) '-'];
	if (allChars == 1) return ' ';

	for (int i = 0; i < 256; i++) if ((char) i != '-' && counts[i] == allChars) return '*';

	if (counts[(unsigned char) 'S'] + 
		counts[(unsigned char) 'T'] + 
		counts[(unsigned char) 'A'] == allChars) 
		return ':';

	if (counts[(unsigned char) 'N'] + 
		counts[(unsigned char) 'E'] + 
		counts[(unsigned char) 'Q'] +
		counts[(unsigned char) 'K'] == allChars) 
		return ':';

	if (counts[(unsigned char) 'N'] + 
		counts[(unsigned char) 'H'] + 
		counts[(unsigned char) 'Q'] +
		counts[(unsigned char) 'K'] == allChars) 
		return ':';

	if (counts[(unsigned char) 'N'] + 
		counts[(unsigned char) 'D'] + 
		counts[(unsigned char) 'E'] +
		counts[(unsigned char) 'Q'] == allChars) 
		return ':';

	if (counts[(unsigned char) 'Q'] + 
		counts[(unsigned char) 'H'] + 
		counts[(unsigned char) 'R'] +
		counts[(unsigned char) 'K'] == allChars) 
		return ':';

	if (counts[(unsigned char) 'M'] + 
		counts[(unsigned char) 'I'] + 
		counts[(unsigned char) 'L'] +
		counts[(unsigned char) 'V'] == allChars) 
		return ':';

	if (counts[(unsigned char) 'M'] + 
		counts[(unsigned char) 'I'] + 
		counts[(unsigned char) 'L'] +
		counts[(unsigned char) 'F'] == allChars) 
		return ':';

	if (counts[(unsigned char) 'H'] + 
		counts[(unsigned char) 'Y'] == allChars) 
		return ':';

	if (counts[(unsigned char) 'F'] + 
		counts[(unsigned char) 'Y'] + 
		counts[(unsigned char) 'W'] == allChars) 
		return ':';

	if (counts[(unsigned char) 'C'] + 
		counts[(unsigned char) 'S'] + 
		counts[(unsigned char) 'A'] == allChars) 
		return '.';

	if (counts[(unsigned char) 'A'] + 
		counts[(unsigned char) 'T'] + 
		counts[(unsigned char) 'V'] == allChars) 
		return '.';

	if (counts[(unsigned char) 'S'] + 
		counts[(unsigned char) 'A'] + 
		counts[(unsigned char) 'G'] == allChars) 
		return '.';

	if (counts[(unsigned char) 'S'] + 
		counts[(unsigned char) 'T'] + 
		counts[(unsigned char) 'N'] + 
		counts[(unsigned char) 'K'] == allChars) 
		return '.';

	if (counts[(unsigned char) 'S'] + 
		counts[(unsigned char) 'T'] + 
		counts[(unsigned char) 'P'] + 
		counts[(unsigned char) 'A'] == allChars) 
		return '.';

	if (counts[(unsigned char) 'S'] + 
		counts[(unsigned char) 'G'] + 
		counts[(unsigned char) 'N'] + 
		counts[(unsigned char) 'D'] == allChars) 
		return '.';

	if (counts[(unsigned char) 'S'] + 
		counts[(unsigned char) 'N'] + 
		counts[(unsigned char) 'D'] + 
		counts[(unsigned char) 'E'] + 
		counts[(unsigned char) 'Q'] + 
		counts[(unsigned char) 'K'] == allChars) 
		return '.';

	if (counts[(unsigned char) 'N'] + 
		counts[(unsigned char) 'D'] + 
		counts[(unsigned char) 'E'] + 
		counts[(unsigned char) 'Q'] + 
		counts[(unsigned char) 'H'] + 
		counts[(unsigned char) 'K'] == allChars) 
		return '.';

	if (counts[(unsigned char) 'N'] + 
		counts[(unsigned char) 'E'] + 
		counts[(unsigned char) 'H'] + 
		counts[(unsigned char) 'Q'] + 
		counts[(unsigned char) 'R'] + 
		counts[(unsigned char) 'K'] == allChars) 
		return '.';

	if (counts[(unsigned char) 'F'] + 
		counts[(unsigned char) 'V'] + 
		counts[(unsigned char) 'L'] + 
		counts[(unsigned char) 'I'] + 
		counts[(unsigned char) 'M'] == allChars) 
		return '.';

	if (counts[(unsigned char) 'H'] + 
		counts[(unsigned char) 'F'] + 
		counts[(unsigned char) 'Y'] == allChars) 
		return '.';

	return ' ';
}

void MultiSequence::WriteALN (std::ostream &outfile, int numColumns)
{
	outfile << "QuickProbs version X.X"  << " multiple sequence alignment" << std::endl;

	int longestComment = 0;
	std::vector<std::vector<char>::iterator> ptrs (count());
	std::vector<int> lengths (count());
	for (int i = 0; i < count(); i++){
		ptrs[i] = GetSequence (i)->getIterator();
		lengths[i] = GetSequence (i)->GetLength();
		longestComment = max (longestComment, (int) GetSequence(i)->getName().length());
	}
	longestComment += 4;

	int writtenChars = 0;    
	bool allDone = false;

	while (!allDone){
		outfile << std::endl;
		allDone = true;

		// loop through all sequences and write them out
		for (int i = 0; i < count(); i++){

			if (writtenChars < lengths[i]){
				outfile << GetSequence(i)->getName();
				for (int j = 0; j < longestComment - (int) GetSequence(i)->getName().length(); j++)
					outfile << ' ';

				for (int j = 0; j < numColumns; j++){
					if (writtenChars + j < lengths[i])
						outfile << ptrs[i][writtenChars + j + 1];
					else
						break;
				}

				outfile << endl;

				if (writtenChars + numColumns < lengths[i]) allDone = false;
			}
		}

		// write annotation line
		for (int j = 0; j < longestComment; j++)
			outfile << ' ';

		for (int j = 0; j < numColumns; j++){
			std::vector<char> column;

			for (int i = 0; i < count(); i++)
				if (writtenChars + j < lengths[i])
					column.push_back (ptrs[i][writtenChars + j + 1]);

			if (column.size() > 0)
				outfile << GetAnnotationChar (column);	
		}

		outfile << endl;
		writtenChars += numColumns;
	}
}

void MultiSequence::SortByHeader () {
	
	// a quick and easy O(n^2) sort
	for (int i = 0; i < (int) sequences.size()-1; i++){
		for (int j = i+1; j < (int) sequences.size(); j++){
			if (sequences[i]->getHeader() > sequences[j]->getHeader())
				std::swap (sequences[i], sequences[j]);
		}
	}
}

void MultiSequence::SortByLabel () {
	
	// a quick and easy O(n^2) sort
	for (int i = 0; i < (int) sequences.size()-1; i++){
		for (int j = i+1; j < (int) sequences.size(); j++){
			if (sequences[i]->GetSortLabel() > sequences[j]->GetSortLabel())
				std::swap(sequences[i], sequences[j]);
		}
	}
}

void MultiSequence::SaveOrdering () {
	
	for (int i = 0; i < (int) sequences.size(); i++)
		sequences[i]->SetSortLabel (i);
}

std::unique_ptr<MultiSequence> MultiSequence::extractSubset(const std::set<int> &indices) const {
	
	std::vector<std::vector<char>::const_iterator> oldPtrs (indices.size());
	std::vector<std::vector<char> *> newPtrs (indices.size());

	assert (indices.size() != 0);

	// grab old data
	int i = 0;
	for (set<int>::const_iterator iter = indices.begin(); iter != indices.end(); ++iter){
		oldPtrs[i++] = GetSequence(*iter)->getIterator();
	}

	// compute new length
	int oldLength = GetSequence (*indices.begin())->GetLength();
	int newLength = 0;
	for (i = 1; i <= oldLength; i++){

		// check to see if there is a gap in every sequence of the set
		bool found = false;
		for (int j = 0; !found && j < (int) indices.size(); j++)
			found = (oldPtrs[j][i] != '-');

		// if not, then this column counts towards the sequence length
		if (found) newLength++;
	}

	// build new alignments
	for (i = 0; i < (int) indices.size(); i++){
		newPtrs[i] = new std::vector<char>(); assert (newPtrs[i]);
		newPtrs[i]->push_back ('@');
	}

	// add all needed columns
	for (i = 1; i <= oldLength; i++){

		// make sure column is not gapped in all sequences in the set
		bool found = false;
		for (int j = 0; !found && j < (int) indices.size(); j++)
			found = (oldPtrs[j][i] != '-');

		// if not, then add it
		if (found){
			for (int j = 0; j < (int) indices.size(); j++)
				newPtrs[j]->push_back (oldPtrs[j][i]);
		}
	}

	// wrap sequences in MultiSequence object
	auto ret = std::unique_ptr<MultiSequence>(new MultiSequence());
	i = 0;
	for (set<int>::const_iterator iter = indices.begin(); iter != indices.end(); ++iter){
		ret->AddSequence (new Sequence (newPtrs[i++], GetSequence (*iter)->getHeader(), newLength,
			GetSequence (*iter)->GetSortLabel(), GetSequence (*iter)->GetLabel()));
	}

	return ret;
}

const ::size_t MultiSequence::calculateHash() const
{
	std::ostringstream oss;
	this->WriteMFA(oss);
	std::hash<std::string> alignmentHash;
	std::string alignmentStr = oss.str();
	return alignmentHash(alignmentStr);
}




