#include <string>
#include <fstream>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include "Sequence.h"

using namespace quickprobs;

Sequence::Sequence(std::vector<char> *data, std::string header, int length, int sequenceLabel, int inputLabel) :
	isValid(data != NULL), header(header), data(data), length (length), sequenceLabel (sequenceLabel), inputLabel (inputLabel) 
{
	assert (data);
	assert ((*data)[0] == '@');
}

Sequence::~Sequence ()
{
	if (data){
		assert (isValid);
		delete data;
		data = NULL;
		isValid = false;
	}
}

Sequence* Sequence::Clone () const 
{
	Sequence *ret = new Sequence();
	assert (ret);

	ret->isValid = isValid;
	ret->header = header;
	ret->data = new std::vector<char>; assert (ret->data);
	*(ret->data) = *data;
	ret->length = length;
	ret->sequenceLabel = sequenceLabel;
	ret->inputLabel = inputLabel;

	return ret;
}

	
Sequence *Sequence::GetRange (int start, int end) const 
{
	Sequence *ret = new Sequence();
	assert (ret);

	assert (start >= 1 && start <= length);
	assert (end >= 1 && end <= length);
	assert (start <= end);

	ret->isValid = isValid;
	ret->header = header;
	ret->data = new std::vector<char>; assert (ret->data);
	ret->data->push_back ('@');
	for (int i = start; i <= end; i++)
		ret->data->push_back ((*data)[i]);
	ret->length = end - start + 1;
	ret->sequenceLabel = sequenceLabel;
	ret->inputLabel = inputLabel;

	return ret;
}

	
Sequence* Sequence::AddGaps (std::vector<char> *alignment, char id) const
{
	Sequence *ret = new Sequence();
	assert (ret);

	ret->isValid = isValid;
	ret->header = header;
	ret->data = new std::vector<char>; assert (ret->data);
	ret->length = (int) alignment->size();
	ret->sequenceLabel = sequenceLabel;
	ret->inputLabel = inputLabel;
	ret->data->push_back ('@');

	std::vector<char>::iterator dataIter = data->begin() + 1;
	for (std::vector<char>::iterator iter = alignment->begin(); iter != alignment->end(); ++iter){
		if (*iter == 'B' || *iter == id){
			ret->data->push_back (*dataIter);
			++dataIter;
		}
		else
			ret->data->push_back ('-');
	}

	return ret;
}

	
std::string Sequence::GetString ()
{
	std::string s = "";
	for (int i = 1; i <= length; i++){
		if ((*data)[i] != '-') s += (*data)[i];
	}
	return s;
}

std::unique_ptr<std::vector<int>> Sequence::getMapping () const 
{
	auto ret = std::unique_ptr<std::vector<int>>(new std::vector<int>(1,0));
	for (int i = 1; i <= length; i++){
		if ((*data)[i] != '-') ret->push_back (i);
	}
	return ret;
}

void Sequence::getMapping(int* ret) const 
{
	ret[0] = 0;
	int i, j;
	for (i = 1, j = 1; i <= length; i++){
		if ((*data)[i] != '-') { 
			ret[j++] = i;
		}
	}
}

	
void Sequence::Highlight (const std::vector<float> &scores, const float cutoff){
	for (int i = 1; i <= length; i++){
		if (scores[i-1] >= cutoff)
			(*data)[i] = toupper ((*data)[i]);
		else
			(*data)[i] = tolower ((*data)[i]);
	}
}

