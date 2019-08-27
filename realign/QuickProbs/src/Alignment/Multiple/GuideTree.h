#pragma once
#include <vector>
#include "Node.h"
#include "Common/StatisticsProvider.h"
#include "Common/Array.h"

//a tree node is a leaf or a node
namespace quickprobs 
{

struct functor_joinBranches;

class GuideTree : public StatisticsProvider
{
	friend struct functor_joinBranches;

public:
	GuideTree(int numSeqs);
	virtual ~GuideTree();	//abstract class

	//get the tree nodes
	const std::vector<Node>& getNodes() const { return nodes; }

	//get the leaf nodes
	const Node* getLeafs() const { return leafs; }
	
	//get the root of the tree
	const Node* getRoot() const { return this->root;}
	Node* getRoot() { return this->root;}

	std::vector<float>& getWeights() { return weights; }

	std::vector<float> calculateSeqsWeights();

	Array<float> calculateBranchDistances();

	Array<float> calculateSubtreeDistances();

	// creates the tree
	virtual void operator()();

	//display the tree
	void displayTree();

protected:
	//system configurations
	int	numSeqs;
	Node* root;
	Node* leafs;
	
	std::vector<float> weights;

	//all the tree nodes
	std::vector<Node> nodes;

	virtual void build() = 0;

	//join two nodes
	void connectNodes(
		Node* parent, int parentIdx, 
		Node* leftChild, float leftDist, 
		Node* rightChild, float rightDist);

};

};