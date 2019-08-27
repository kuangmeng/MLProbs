#include "AminoAcidMatrices.hpp"

// integer specialisations
template class Blosum45<int>;
template class Blosum50<int>;
template class Blosum62<int>;
template class Blosum80<int>;

template class Pam30<int>;
template class Pam70<int>;
template class Pam120<int>;
template class Pam250<int>;

//template class Gonnet160<int>;
//template class Vtml200<int>;

// float specialisations
template class Blosum45<float>;
template class Blosum50<float>;
template class Blosum62<float>;
template class Blosum80<float>;

template class Pam30<float>;
template class Pam70<float>;
template class Pam120<float>;
template class Pam250<float>;

template class Gonnet160<float>;
template class Vtml200<float>;


// double specialisations
template class Blosum45<double>;
template class Blosum50<double>;
template class Blosum62<double>;
template class Blosum80<double>;

template class Pam30<double>;
template class Pam70<double>;
template class Pam120<double>;
template class Pam250<double>;

template class Gonnet160<double>;
template class Vtml200<double>;

template class Miqs<double>;
template class MiqsFP<double>;