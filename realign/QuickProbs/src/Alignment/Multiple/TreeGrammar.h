#pragma once
#include <boost/spirit/home/classic.hpp>
#include <boost/spirit/home/classic/phoenix.hpp>

#include "GuideTree.h"
#include "Alignment/DataStructures/Sequence.h"
#include "Alignment/DataStructures/MultiSequence.h"

namespace bs = boost::spirit::classic;
namespace ph = phoenix;
namespace phs = phoenix;

namespace quickprobs
{

// forward declarations


struct Branch
{
	int id;
	double length;

	Branch() : id(0), length(0) {}
	Branch(int id, double length) : id(id), length(length) {}
};

struct TreeWrapper {
    GuideTree& tree;
	std::map<std::string, int> sequencesToIds;
	int currentInternalId;

	TreeWrapper(GuideTree& tree, MultiSequence& sequences) : tree(tree) {
        // fill in mappings
		for (int i = 0; i < sequences.count(); ++i) {
			auto seq = sequences.GetSequence(i);
			sequencesToIds[seq->getName()] = i;
		}

		currentInternalId = sequences.count();
	}
};

struct functor_setLength
{
	template <typename T1, typename T2>
	struct result {	typedef void type; };

	template <typename T1, typename T2>
	void operator()(T1& branch, T2 length) const { branch.length = length; }
};


struct functor_sequenceToId
{
	TreeWrapper& wrapper;
	functor_sequenceToId(TreeWrapper &wrapper) : wrapper(wrapper) {}

	template <typename T>
	struct result {	typedef Branch type; };

	template <typename T>
	Branch operator()(T name) const
	{
		int id =  wrapper.sequencesToIds[name];
		//std::cout << "sequences[" << id << "] = " << name << std::endl;
		Branch out(id, 0.0);
		return out;
	}
};


struct functor_joinBranches
{
	TreeWrapper& wrapper;
	functor_joinBranches(TreeWrapper &wrapper) : wrapper(wrapper) {}

	template <typename T1, typename T2>
	struct result { typedef Branch type; };

	template <typename T1, typename T2>
	Branch operator()(T1 left, T2 right) const
	{
		float leftLength = left.length;
		float rightLength = right.length;

		Node* leftChild = &wrapper.tree.nodes[left.id];
		Node* rightChild = &wrapper.tree.nodes[right.id];
		Node* parent = &wrapper.tree.nodes[wrapper.currentInternalId];

		wrapper.tree.connectNodes(parent, wrapper.currentInternalId, leftChild, left.length, rightChild, right.length);

	//	std::cout << "getNodes()[" << wrapper.currentInternalId << "] = join(" <<
	//		leftChild->idx << ":" << leftChild->dist << "," << rightChild->idx << ":" << rightChild->dist << ")" << std::endl;

		Branch out(wrapper.currentInternalId++, 0.0);
		return out;
	}
};

template <typename T>
struct SingleClosure : bs::closure<SingleClosure<T>, T>
{
    typename bs::closure<SingleClosure<T>, T>::member1 val;
};

struct TreeGrammar : public bs::grammar<TreeGrammar>
{
	TreeWrapper wrapper;

	ph::function<functor_sequenceToId> sequenceToId;
	ph::function<functor_joinBranches> joinBranches;
	ph::function<functor_setLength> setLength;

	TreeGrammar(GuideTree& tree, MultiSequence& sequences) :
        wrapper(tree, sequences),
		sequenceToId(functor_sequenceToId(wrapper)),
		joinBranches(functor_joinBranches(wrapper))
	{
	}


	template <typename ScannerT>
	struct definition
	{
		bs::rule<ScannerT, SingleClosure<std::string>::context_t> seqname;
		bs::rule<ScannerT, SingleClosure<double>::context_t> length;

		bs::rule<ScannerT, SingleClosure<Branch>::context_t> node;
		bs::rule<ScannerT, SingleClosure<Branch>::context_t> branch;
		bs::rule<ScannerT, SingleClosure<Branch>::context_t> branchSet;
		bs::rule<ScannerT> head;

		definition(TreeGrammar const& self)
		{
			// grammar head
			head = node >> ";";

			// non-terminal tokens
			branchSet = branch[branchSet.val = phs::arg1] >> +("," >> branch[branchSet.val = self.joinBranches(branchSet.val, phs::arg1)]);

			branch = node[branch.val = phs::arg1] >> ":" >> length[self.setLength(branch.val, phs::arg1)];

			node = seqname[node.val = self.sequenceToId(phs::arg1)]	|
				("(" >> branchSet[node.val = phs::arg1] >> ")") ;

			// terminal tokens (sequence names or lengths)
			seqname = (+(bs::alnum_p | '_'))[seqname.val = ph::construct_<std::string>(phs::arg1, phs::arg2)];
			length = bs::real_p[length.val = phs::arg1];

		};

		bs::rule<ScannerT> const& start() const { return head; }
	};
};

};
