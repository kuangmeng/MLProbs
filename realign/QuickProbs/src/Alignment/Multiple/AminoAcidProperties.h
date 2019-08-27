#pragma once
#include <vector>

class AminoAcidProperties
{
public:
	static const unsigned int SMALL;
	static const unsigned int TINY;
	static const unsigned int ALIPHATIC;
	static const unsigned int AROMATIC;
	static const unsigned int HYDROPHOBIC;
	static const unsigned int POLAR;
	static const unsigned int POSITIVE;
	static const unsigned int NEGATIVE;
	static const unsigned int CHARGED;
	static const unsigned int PROLINE;
	
	AminoAcidProperties();
	unsigned int get(char a) { return props[a]; }
protected:
	std::vector<unsigned int> props;
};