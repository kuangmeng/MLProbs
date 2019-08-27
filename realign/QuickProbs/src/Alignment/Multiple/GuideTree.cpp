#include <iostream>
#include "Common/Timer.h"
#include "GuideTree.h"

#undef max

using namespace quickprobs;


GuideTree::GuideTree(int numSeqs) : numSeqs(numSeqs)
{
	int i;
	Node* node;
	
	//system configuration
	nodes.resize(numSeqs * 2 + 1);

	//initialize all the tree nodes
	leafs = nodes.data();
	
	for (i = 0; i < nodes.size(); i ++){
		node = &nodes[i];
		node->left = nullptr;
		node->right = nullptr;
		node->parent = nullptr; 
		node->leftIdx = -1; 
		node->rightIdx = -1; 
		node->parentIdx = -1;
		node->idx = -1; 
		node->dist = 0;
		node->type = INTERNAL_NODE;		//set to be NODE, by default
		node->order = 0;
		node->depth = 0;
	}
	
	//initialize the leaf nodes
	for (i = 0; i < numSeqs; i++){
		node = &leafs[i];
		node->idx = i;
		node->type = LEAF_NODE;
	}
}

GuideTree::~GuideTree()
{
}

void GuideTree::operator()()
{
	TIMER_CREATE(timer);
	TIMER_START(timer);

	this->build();
	this->weights = calculateSeqsWeights();
	
	TIMER_STOP(timer);
	int height = std::max_element(leafs, leafs + numSeqs, [](const Node& n1, const Node &n2)->int { 
		return n1.depth < n2.depth; 
	})->depth; 
}

/****************************************************
	create the evolutionary relationship
****************************************************/
void GuideTree::connectNodes(
	Node* parent, int parentIdx, 
	Node* leftChild, float leftDist, 
	Node* rightChild, float rightDist)
{
	//save the parents index for each child
	leftChild->parent = parent;
	leftChild->parentIdx = parentIdx;
	rightChild->parent  = parent;
	rightChild->parentIdx = parentIdx;

	//save the branch lengths (i.e. distance) from each child to its parent
	leftChild->dist = leftDist;
	rightChild->dist = rightDist;

	//save the indices of itself and its children for this new tree node
	parent->idx = parentIdx;
	parent->left = leftChild;
	parent->leftIdx = leftChild->idx;
	parent->right = rightChild;
	parent->rightIdx = rightChild->idx;

}

/*********************************
	display the tree
*********************************/
void GuideTree::displayTree()
{	
	fprintf(stderr, "**************DISPLAY TREE*********************\n");
	for(int i = 0; i < nodes.size(); i++){
		Node* node = &nodes[i];

		fprintf(stderr, "%d(%p): left(%p) %d, right(%p) %d, parent(%p) %d, dist %f\n",
			(node == &nodes[node->idx]) ? node->idx: -2,
			node,
			node->left,	
			(!node->left || node->left == &nodes[node->leftIdx]) ? node->leftIdx : -2,
			node->right,
			(!node->right || node->right == &nodes[node->rightIdx]) ? node->rightIdx : -2,
			node->parent,
			(!node->parent || node->parent == &nodes[node->parentIdx]) ? node->parentIdx : -2, node->dist);
	}
	fprintf(stderr, "*******************************************\n");
}

/*********************************
	compute the sequence weights
*********************************/
std::vector<float> GuideTree::calculateSeqsWeights()
{
	std::vector<float> weights(numSeqs);
	
	//compute the order of each node, which represents the number of leaf nodes in the sub tree rooting from it.
	std::for_each(leafs, leafs + numSeqs, [](Node& leaf)->void {
		Node* curr = &leaf;
		while(curr != nullptr) {
			++curr->order;
			++leaf.depth;
			curr = curr->parent;
		} 
	});
	
	//compute the weight of each sequence, which corresponds to a leaf node
	std::transform(leafs, leafs + numSeqs, weights.begin(), [](const Node& leaf)->float {
		const Node* curr = &leaf;
		float w = 0;
		while (curr->parent != nullptr) {
			w += curr->dist / curr->order;
			curr = curr->parent;
		}
		//return (float)(int)(w * 100);
		return w;
	}); 

	// calculate sum of weights
	float wsum = std::accumulate(weights.begin(), weights.end(), 0.0f);
	if (wsum == 0) {
		//in this case, every sequence is assumed to have an identical weight
		std::fill(weights.begin(), weights.end(), 1.0f); 
		wsum = numSeqs;
	} 
	
	// normalize weights
	std::transform(weights.begin(), weights.end(), weights.begin(), [wsum](float w)->float {
		return w / wsum;
	});

	return weights;
}


Array<float> quickprobs::GuideTree::calculateBranchDistances()
{
	Array<float> treeDists(numSeqs);
	
	std::vector<std::string> paths(numSeqs);
	for (auto &p : paths) { p.reserve(numSeqs); }

	for (int i = 0; i < numSeqs; ++i) {
		Node *n = &leafs[i];

		while (n->parent != nullptr) {
			if (n == n->parent->left) { paths[i].push_back('l'); } 
			else { paths[i].push_back('r'); }
			n = n->parent;
		}
	}

	for (int i = 0; i < numSeqs; ++i) {
		for (int j = i + 1; j < numSeqs; ++j) {
			auto& p1 = (paths[i].size() > paths[j].size()) ? paths[i] : paths[j];
			auto& p2 = (paths[i].size() > paths[j].size()) ? paths[j] : paths[i]; 
			// check common suffixes
			auto it = std::mismatch(p1.rbegin(), p1.rend(), p2.rbegin());
			size_t d1 = p1.rend() - it.first;
			size_t d2 = p2.rend() - it.second;
			treeDists[i][j] = treeDists[j][i] = d1 + d2;
		}
	}

	return treeDists;
}

Array<float> quickprobs::GuideTree::calculateSubtreeDistances()
{
	Array<float> treeDists(numSeqs);
	std::vector<std::vector<int>> paths(numSeqs);
	
	for (auto &p : paths) { p.reserve(numSeqs); }

	for (int i = 0; i < numSeqs; ++i) {
		Node *n = &leafs[i];
		
		do {
			paths[i].push_back(n->idx);
			n = n->parent;
		} while (n != nullptr);
	}

	for (int i = 0; i < numSeqs; ++i) {
		for (int j = i + 1; j < numSeqs; ++j) {
			auto& p1 = (paths[i].size() > paths[j].size()) ? paths[i] : paths[j];
			auto& p2 = (paths[i].size() > paths[j].size()) ? paths[j] : paths[i]; 
			// check common suffixes
			auto it = std::mismatch(p1.rbegin(), p1.rend(), p2.rbegin());
			
			// iterate through common part
			auto& id1 = *(it.first);
			auto& id2 = *(it.second);
			size_t dist = nodes[id1].order + nodes[id2].order;
			treeDists[i][j] = treeDists[j][i] = dist;
		}
	}

	return treeDists;
}



