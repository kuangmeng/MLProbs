#include "NucleotideMatrices.hpp"

// integer specialisations
template class NaiveNucleotide<int>;
template class Hoxd55<int>;
template class Hoxd70<int>;

// double specialisations
template class NaiveNucleotide<double>;
template class Hoxd55<double>;
template class Hoxd70<double>;
