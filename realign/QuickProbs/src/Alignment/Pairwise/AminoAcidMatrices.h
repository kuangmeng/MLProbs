#pragma once
#include "Scoring.h"

// It is assumed that at the gap opening only gi penalty is given.
// SSEARCH assumes that at gap opening both gi and ge are given. This is why
// we add ge to gi to obtain the same values.

template <class T>
class Blosum45 : public ProteinScoring<T>
{
public:
	static const T SCORES[];
	Blosum45(T ge = -2, T gi = -15) : ProteinScoring<T>(SCORES, ge, gi){}
};

template <class T>
class Blosum50 : public ProteinScoring<T>
{
public:
	static const T SCORES[];
	Blosum50(T ge = -2, T gi = -10) : ProteinScoring<T>(SCORES, ge, gi){}
};

template <class T>
class Blosum62 : public ProteinScoring<T>
{
public:
	static const T SCORES[];
	Blosum62(T ge = -1, T gi = -7) : ProteinScoring<T>(SCORES, ge, gi){}
};

template <class T>
class Blosum80 : public ProteinScoring<T>
{
public:
	static const T SCORES[];
	Blosum80(T ge = -4, T gi = -16) : ProteinScoring<T>(SCORES, ge, gi){}
};

template <class T>
class Pam30 : public ProteinScoring<T>
{
public:
	static const T SCORES[];
	Pam30(T ge = -1, T gi = -9) : ProteinScoring<T>(SCORES, ge, gi){}
};

template <class T>
class Pam70 : public ProteinScoring<T>
{
public:
	static const T SCORES[];
	Pam70(T ge = -1, T gi = -10) : ProteinScoring<T>(SCORES, ge, gi){}
};

template <class T>
class Pam120 : public ProteinScoring<T>
{
public:
	static const T SCORES[];
	Pam120(T ge = -4, T gi = -16) : ProteinScoring<T>(SCORES, ge, gi){}
};

template <class T>
class Pam250 : public ProteinScoring<T>
{
public:
	static const T SCORES[];
	Pam250(T ge = -2, T gi = -2) : ProteinScoring<T>(SCORES, ge, gi){}
};

template <class T>
class Gonnet160 : public ProteinScoring<T>
{
public:
	static const T SCORES[];
	Gonnet160(T ge = -1.0, T gi = -22.0) : ProteinScoring<T>(SCORES, ge, gi) {}
};

template <class T>
class Vtml200 : public ProteinScoring<T>
{
public:
	static const T SCORES[];
	Vtml200(T ge = -1.5, T gi = -22.15) : ProteinScoring<T>(SCORES, ge, gi) {}
};

template <class T>
class Miqs : public ProteinScoring<T>
{
public:
	static const T SCORES[];
	Miqs(T ge = -1.5, T gi = -22.15) : ProteinScoring<T>(SCORES, ge, gi) {}
};

template <class T>
class MiqsFP : public ProteinScoring<T>
{
public:
	static const T SCORES[];
	MiqsFP(T ge = -1.5, T gi = -22.15) : ProteinScoring<T>(SCORES, ge, gi) {}
};