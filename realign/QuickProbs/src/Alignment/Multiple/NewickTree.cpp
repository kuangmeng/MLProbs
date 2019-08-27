#include <stdexcept>
#include <boost/spirit/include/classic.hpp>

#include "DataStructures/MultiSequence.h"
#include "NewickTree.h"
#include "TreeGrammar.h"

namespace bs = boost::spirit::classic;
using namespace quickprobs;

NewickTree::NewickTree(quickprobs::MultiSequence& sequences, std::string description) 
	: GuideTree(sequences.count()), sequences(sequences), description(description)
{	
}

void NewickTree::build()
{
	if (description.length() == 0) {
		throw std::runtime_error("NewickTree::build(): 'description' member not set.");
	}

	TreeGrammar grammar(*this, sequences);
	
	auto info = bs::parse(description.c_str(), grammar);
	int rootId = grammar.wrapper.currentInternalId - 1;
	this->root = &nodes[rootId];

	if (info.full == false) {
		throw std::runtime_error("NewickTree::build(): invalid tree format.");
	}
}

