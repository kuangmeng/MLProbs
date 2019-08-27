#include "AminoAcidProperties.h"


const unsigned int AminoAcidProperties::SMALL = 1;
const unsigned int AminoAcidProperties::TINY = 2;
const unsigned int AminoAcidProperties::ALIPHATIC = 4;
const unsigned int AminoAcidProperties::AROMATIC = 8;
const unsigned int AminoAcidProperties::HYDROPHOBIC = 16;
const unsigned int AminoAcidProperties::POLAR = 32;
const unsigned int AminoAcidProperties::POSITIVE = 64;
const unsigned int AminoAcidProperties::NEGATIVE = 128;
const unsigned int AminoAcidProperties::CHARGED = 256;
const unsigned int AminoAcidProperties::PROLINE = 512;

AminoAcidProperties::AminoAcidProperties()
{
	props.resize(256, 0xffffffff);
	
	props['A'] = TINY | SMALL | HYDROPHOBIC;
	props['C'] = TINY | SMALL | HYDROPHOBIC | POLAR;
	props['D'] = NEGATIVE | CHARGED | POLAR | SMALL;
	props['E'] = NEGATIVE | CHARGED | POLAR;
	props['F'] = AROMATIC | HYDROPHOBIC;
	props['G'] = TINY | SMALL | HYDROPHOBIC;
	props['H'] = AROMATIC | POSITIVE | CHARGED | POLAR | HYDROPHOBIC;
	props['I'] = ALIPHATIC | HYDROPHOBIC; 
	props['K'] = POSITIVE | CHARGED | POLAR | HYDROPHOBIC; 
	props['L'] = ALIPHATIC | HYDROPHOBIC; 
	props['M'] = HYDROPHOBIC; 
	props['N'] = SMALL | POLAR; 
	props['P'] = PROLINE | SMALL; 
	props['Q'] = POLAR;
	props['R'] = POSITIVE | CHARGED | POLAR; 
	props['S'] = TINY | SMALL | POLAR; 
	props['T'] = POLAR | HYDROPHOBIC | SMALL;
	props['V'] = SMALL | ALIPHATIC | HYDROPHOBIC;
	props['W'] = AROMATIC | POLAR | HYDROPHOBIC; 
	props['Y'] = AROMATIC | POLAR | HYDROPHOBIC;
}
