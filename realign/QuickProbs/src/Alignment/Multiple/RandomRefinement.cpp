#include "RandomRefinement.h"

using namespace quickprobs;


void quickprobs::RandomRefinement::split(
	const GuideTree& tree,
	const MultiSequence& alignment, 
	std::set<int>& groupOne, 
	std::set<int>& groupTwo)
{
	int numSeqs = alignment.count();

	int index = genRandom(numSeqs, this->currentIter, true);
	// create two separate groups
	for (int i = 0; i < numSeqs; ++i) {
		index = genRandom(numSeqs, this->currentIter, false);
		if (index % 2) {
			groupOne.insert(i);
		} else {
			groupTwo.insert(i);
		}
	}
}

int quickprobs::RandomRefinement::genRandom(int m, int seed, bool init)
{
	static const int a= 5, b = 3, n = 7;
	static int rand0;
	if(init == true){
		rand0 = seed;
	}
	m *= 19;
	int rand1;
	for(int i = 0; i < n; i++){
		rand1 = (a * rand0 + b) % m;
		rand0 = rand1;
	}
	return rand1;
}


