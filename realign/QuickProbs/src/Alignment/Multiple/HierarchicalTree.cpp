#include "HierarchicalTree.h"

#include <map>
#include <list>
#include <algorithm>

using namespace quickprobs;
using namespace std;

HierarchicalTree::HierarchicalTree(std::vector<std::vector<float>>& distances) :
	GuideTree(distances.size()), distances(distances)
{

}

void HierarchicalTree::build()
{
	std::vector<Node*> bestMerge(numSeqs);

	// generate all edges
	int did = 0;
	for (int i = 0; i < numSeqs; ++i) {
		int min_id = 0;
		float min_distance = std::numeric_limits<float>::max();
		for (int j = 0; j < numSeqs; ++j) {
			min_id = (distances[i][j] < min_id && i != j) ? distances[i][j] : min_distance;
		}
		bestMerge[i] = &leafs[min_id];
	}

	for (int step = 0; step < numSeqs - 1; ++step) {
		auto minN = *std::min_element(bestMerge.begin(), bestMerge.end());
	}
}

