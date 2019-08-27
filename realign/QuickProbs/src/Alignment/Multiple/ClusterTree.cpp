#include <stdexcept>
#include "ClusterTree.h"

#include "Common/Array.h"

using namespace quickprobs;

ClusterTree::ClusterTree(Array<float>& distances) 
	: GuideTree(distances.size()), distances(distances)
{
}

ClusterTree::~ClusterTree()
{
}

void ClusterTree::build()
{
	int i;
	TempNode *headValidNodes;
	TempNode* miniPtr, *minjPtr, *ivalid, *jvalid;
	int mini, minj;

	std::vector<TempNode>validNodes(numSeqs + 1);
	std::vector<float>joins(numSeqs + 1);
	std::vector<unsigned int>clusterLeafs(nodes.size() + 1);

	//initialize cluster size 
	for(i = 0; i < numSeqs; i++){
		clusterLeafs[i] = 1;
	}

	headValidNodes = &validNodes[0];
	headValidNodes->next = &validNodes[1];
	headValidNodes->n = -1;
	headValidNodes->node = -1;
	headValidNodes->prev = NULL;

	//build an initial link list
	TempNode* curr = &validNodes[1];
	TempNode* prev = headValidNodes;
	TempNode* next = &validNodes[2];
	for(i = 0; i < numSeqs; i++){
		curr->n = i;
		curr->node = i;
		curr->prev = prev;
		curr->next = next;
		prev = curr;
		curr = next;
		next++;
	}
	prev->next = NULL;

	//to generate the cluster tree
	int nodeIdx;	//the index of an internal node
	int firstNode = numSeqs;	//the index of the first internal node
	int lastNode = firstNode + numSeqs - 1;	//the index of the last internal node

	for(nodeIdx = firstNode; nodeIdx < lastNode; nodeIdx++){
		//find closest pair of clusters
		float minDist = 2.0f;
		miniPtr = headValidNodes;
		minjPtr = headValidNodes;

		for(ivalid = headValidNodes->next; ivalid != NULL; ivalid = ivalid->next){
			mini = ivalid->n;
			for(jvalid = headValidNodes->next; jvalid != NULL && jvalid->n < mini; jvalid = jvalid->next){
				minj = jvalid->n;
				float dist = distances[mini][minj];
				if(dist < 0){
					throw std::runtime_error("ERROR: It is impossible to have distance value less than zero");
					dist = 0;
				}
				if(dist < minDist){
					minDist = dist;
					miniPtr = ivalid;
					minjPtr = jvalid;
				}
				//printf("dist %g mini %d minj %d\n", dist, ivalid->node, jvalid->node);
			}
		}
		//printf("**** mini %d minj %d minDist %g *****\n", miniPtr->node, minjPtr->node, minDist);
		//check the validity of miniPtr and minjPtr;
		if(miniPtr == headValidNodes || minjPtr == headValidNodes){
			throw std::runtime_error("OOPS: Error occurred while constructing the cluster tree\n");
		}
		//computing branch length and join the two nodes
		float branchLength = minDist * 0.5f;
		this->connectNodes(&nodes[nodeIdx], nodeIdx, &nodes[miniPtr->node], branchLength,
			&nodes[minjPtr->node], branchLength);
		clusterLeafs[nodeIdx] = clusterLeafs[miniPtr->node] + clusterLeafs[minjPtr->node];

		//remove the valid node minjPtr from the list
		minjPtr->prev->next = minjPtr->next;
		if(minjPtr->next != NULL){
			minjPtr->next->prev = minjPtr->prev;
		}
		minjPtr->prev = minjPtr->next = NULL;

		//compute the distance of each remaining valid node to the new node
		for(ivalid = headValidNodes->next; ivalid != NULL; ivalid = ivalid->next){
			int idx = ivalid->n;

			float idist = distances[miniPtr->n][idx];
			float jdist = distances[minjPtr->n][idx];

			unsigned int isize = clusterLeafs[miniPtr->node];
			unsigned int jsize = clusterLeafs[minjPtr->node];
			joins[idx] = (idist * isize + jdist * jsize) / (isize + jsize);
		}
		//update the distance to the new node
		miniPtr->node = nodeIdx;
		mini = miniPtr->n;
		for(jvalid = headValidNodes->next; jvalid != NULL; jvalid = jvalid->next){
			minj = jvalid->n;

			float dist = joins[minj];
			distances[mini][minj] = dist;
			distances[minj][mini] = dist;
		}
	}
	//add a pseudo root to this unrooted NJ tree
	this->root = &nodes[lastNode - 1];
}
