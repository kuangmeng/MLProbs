#pragma once
#include "Scoring.h"

template <class T>
class NaiveNucleotide : public DnaScoring<T>
{
public:
	static const T SCORES[];
	NaiveNucleotide(T ge = -1.0, T gi = -7.0) : DnaScoring<T>(SCORES, ge, gi){}
};

template <class T>
class Hoxd55 : public DnaScoring<T>
{
public:
	static const T SCORES[];
	Hoxd55(T ge = -30.0, T gi = -400.0) : DnaScoring<T>(SCORES, ge, gi){}
};

template <class T>
class Hoxd70 : public DnaScoring<T>
{
public:
	static const T SCORES[];
	Hoxd70(T ge = -30.0, T gi = -400.0) : DnaScoring<T>(SCORES, ge, gi){}
};