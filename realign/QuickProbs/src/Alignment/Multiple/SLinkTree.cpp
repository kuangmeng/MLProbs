#include "SLinkTree.h"

#include <vector>
#include <algorithm>
#include <functional>

using namespace std;
using namespace quickprobs;

SLinkTree::SLinkTree(Array<float>& distances) 
	: GuideTree(distances.size()), distances(distances)
{
}

void SLinkTree::build()
{
	int i, j;
	int next;
	
	vector<int> pi(numSeqs);
	vector<float> lambda(numSeqs);
	vector<float> sim_vector(numSeqs);

	for(i = 0; i < numSeqs; ++i) {
		pi[i] = i;
		lambda[i] = -std::numeric_limits<float>::max();

/*#ifdef SHOW_PROGRESS
		if(n_seq > 1000 && i % 100 == 0)
		{
			cout << i << " / " << n_seq << "\r";
			fflush(stdout);
		}
#endif*/

		std::transform(distances[i], distances[i] + numSeqs, sim_vector.begin(), [](float x)->float{
			return 1 - x;
		});

		for(j = 0; j < i; ++j)
		{
			next = pi[j];
			if(lambda[j] > sim_vector[j])
				sim_vector[next] = max(sim_vector[next], sim_vector[j]);
			else
			{
				sim_vector[next] = max(lambda[j], sim_vector[next]);
				pi[j] = i;
				lambda[j] = sim_vector[j];
			}
		}

		for(j = 0; j < i; ++j)	
		{
			next = pi[j];
			if(lambda[next] >= lambda[j])
				pi[j] = i;
		}
	}

	vector<int> elements(numSeqs-1);
	for(i = 0; i < numSeqs-1; ++i)
		elements[i] = i;

	vector<int> index(numSeqs);
	for(i = 0; i < numSeqs; ++i) {
		index[i] = i;
	}

	stable_sort(elements.begin(), elements.end(), [&](int x, int y){
		return lambda[x] > lambda[y];
	});

	stable_sort(lambda.begin(), lambda.end(), std::greater<int>());

/*	std::vector<std::pair<int,int>> guide_tree;
	for(size_t i = 0; i < numSeqs; ++i) {
		guide_tree.push_back(make_pair(-1, -1));
	}

	for (i = 0; i < numSeqs-1; ++i) {
		j = elements[i];
		next = pi[j];
		guide_tree.push_back(make_pair(index[j], index[next]));
		index[next] = numSeqs + i;
	}

	int a = 0;*/

	
	for (i = 0; i < numSeqs-1; ++i) {
		j = elements[i];
		next = pi[j];

		int leftIndex = index[j];
		int rightIndex = index[next];

		Node * leftChild = &nodes[leftIndex];
		Node * rightChild = &nodes[rightIndex];
		Node * parent = &nodes[i + numSeqs];

		float branchLength = (1 - lambda[i]) / 2.0f;
		this->connectNodes(parent, i + numSeqs, leftChild, branchLength ,rightChild, branchLength);
		auto p = make_pair(index[j], index[next]);

		index[next] = numSeqs + i;
	} 

	this->root = &nodes[2 * numSeqs - 2];
}
