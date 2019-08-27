#include "TreeRefinement.h"

#include "Common/deterministic_random.h"
#include <vector>
#include <functional>

using namespace std;
using namespace quickprobs;

void TreeRefinement::split(
	const GuideTree& tree,	
	const MultiSequence& alignment,
	std::set<int>& groupOne,
	std::set<int>& groupTwo)
{
	int numSeqs = alignment.count();

	vector<bool> inGroup1(numSeqs);
	for (int i = 0; i < numSeqs; i++) {
		inGroup1[i] = false;
	}

	int rootId = tree.getRoot()->idx;
	det_uniform_int_distribution<int> dist(0, rootId-1); // exclude root
	int nodeId = dist(engine);

	Node node = tree.getNodes()[nodeId];

	std::function<void(const Node&)> fun = [&](const Node& node) {
		if (node.type == NodeType::LEAF_NODE) {
			inGroup1[node.idx] = true;
		}else {
			fun(*node.left);
			fun(*node.right);
		}
	};

	fun(node);

	// create two separate groups
	for (int i = 0; i < numSeqs; i++) {
		if (inGroup1[i]) {
			groupOne.insert(i);
		} else {
			groupTwo.insert(i);
		}
	}
}