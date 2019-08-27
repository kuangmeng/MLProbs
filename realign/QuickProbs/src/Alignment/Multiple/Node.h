#pragma once

//a tree node is a leaf or a node
namespace quickprobs 
{

	enum NodeType {
		EMPTY_NODE,
		INTERNAL_NODE,
		LEAF_NODE
	};

	struct TempNode
	{
		TempNode*	prev;
		TempNode*	next;
		int			n;				//the index in the distance matrix			
		int			node;			//the index in the tree node entries
	};

	struct Node
	{
		struct Node *left;			//the pointer to its left child
		struct Node *right;			//the pointer to its right child
		struct Node *parent;		//the pointer to its parent
		int leftIdx;					//the index of the left child
		int rightIdx;					//the index of the right child
		int parentIdx;					//the index of its parent
		int idx;						//the index of itself
		float dist;						//the distance to its parent
		NodeType type;					//whether it is a leaf node or not
		int order;						//the number of generations dating back to its ancestor
		int	depth;						//the depth of the node
	};
}