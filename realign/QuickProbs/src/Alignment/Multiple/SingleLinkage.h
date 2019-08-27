#include "GuideTree.h"

namespace quickprobs
{
class SingleLinkage : public GuideTree {
public:
		SingleLinkage(std::vector<std::vector<float>>& distances) : 
			GuideTree(distances.size()), distances(distances) {}
		
		// creates the tree
		virtual void build();

protected:
		std::vector<std::vector<float>>& distances;
};
}