#pragma once
#include<stdexcept>

class TreeKind {
public:
	enum Enumeration {UPGMA, SingleLinkage, Chained};
	TreeKind() {}
	TreeKind(std::string s) {
		if (s == "upgma") {  value = UPGMA; }
		else if (s == "slink") { value = SingleLinkage; }
		else if (s == "chained") { value = Chained; }
		else { throw std::runtime_error("Illegal TreeKind!"); } 
	}

	std::string toString() {
		switch (value) {
			case UPGMA: return "upgma";
			case SingleLinkage: return "slink";
			case Chained: return "chained";
		}
		return "invalid";
	}

	bool operator==(TreeKind::Enumeration ref) { return value == ref; }
	TreeKind& operator=(TreeKind::Enumeration ref) { value = ref; return *this; }

protected:
	Enumeration value;
};


class SelectivityFunction
{
public:
	enum Enumeration {Sum, Min, Max, Avg};

	SelectivityFunction() {}
	SelectivityFunction(std::string s) {
		if (s == "sum") {  value = Sum; }
		else if (s == "min") { value = Min; }
		else if (s == "max") { value = Max; }
		else if (s == "avg") { value = Avg; }
		else { throw std::runtime_error("Illegal SelectivityTransformation!"); } 
	}

	std::string toString() {
		switch (value) {
			case Sum: return "sum";
			case Min: return "min";
			case Max: return "max";
			case Avg: return "avg";
		}
		return "invalid";
	}

	bool operator==(SelectivityFunction::Enumeration ref) { return value == ref; }
	SelectivityFunction& operator=(SelectivityFunction::Enumeration ref) { value = ref; return *this; }

protected:
	Enumeration value;
};

class SelectivityNormalization {
public:
	enum Enumeration {No, Stochastic, RankedStochastic, RankedRowStochastic};

	SelectivityNormalization() {}
	SelectivityNormalization(std::string s) {
		if (s == "no") { value = No; }
		else if (s == "stochastic") { value = Stochastic; }
		else if (s == "ranked-stochastic") { value = RankedStochastic; }
		else if (s == "ranked-row-stochastic") { value = RankedRowStochastic; }
		else { throw std::runtime_error("Illegal SelectivityNormalization!"); } 
	}

	std::string toString() {
		switch (value) {
			case SelectivityNormalization::No: return "no";
			case SelectivityNormalization::Stochastic: return "stochastic";
			case SelectivityNormalization::RankedStochastic: return "ranked-stochastic";
			case SelectivityNormalization::RankedRowStochastic: return "ranked-row-stochastic";
		}
		return "invalid";
	}

	bool operator==(SelectivityNormalization::Enumeration ref) { return value == ref; }
	SelectivityNormalization& operator=(SelectivityNormalization::Enumeration ref) { value = ref; return *this; }

protected:
	Enumeration value;
};


class SelectivityMode {
public:
	enum Enumeration {Subtree, Seed, Similarity};

	SelectivityMode() {}
	SelectivityMode(std::string s) {
		if (s == "subtree") { value = SelectivityMode::Subtree; }
		else if (s == "seed") { value = SelectivityMode::Seed; }
		else if (s == "similarity") { value = SelectivityMode::Similarity; }
		else { throw std::runtime_error("Illegal SelectivityMode!"); }
	}

	std::string toString() {
	switch (value) {
			case SelectivityMode::Subtree : return "subtree";
			case SelectivityMode::Seed: return "seed";
			case SelectivityMode::Similarity: return "similarity";
		} 
		return "invalid";
	}

	bool operator==(SelectivityMode::Enumeration ref) { return value == ref; }
	SelectivityMode& operator=(SelectivityMode::Enumeration ref) { value = ref; return *this; }

protected:
	Enumeration value;
};




class SelectivityFilter {
public:
	enum Enumeration {Deterministic, TriangleLowpass, TriangleMidpass, TriangleHighpass, HomographLowpass};

	SelectivityFilter() {}
	SelectivityFilter(std::string s) {
		if (s == "deterministic") { value = SelectivityFilter::Deterministic; }
		else if (s == "triangle-lowpass") { value = SelectivityFilter::TriangleLowpass; }
		else if (s == "triangle-midpass") { value = SelectivityFilter::TriangleMidpass; }
		else if (s == "triangle-highpass") { value = SelectivityFilter::TriangleHighpass; }
		else if (s == "homograph-lowpass") { value = SelectivityFilter::HomographLowpass; }
		else { throw std::runtime_error("Illegal SelectivityDistribution!"); }
	}

	std::string toString() {
	switch (value) {
			case SelectivityFilter::Deterministic: return "deterministic";
			case SelectivityFilter::TriangleLowpass : return "triangle-lowpass";
			case SelectivityFilter::TriangleMidpass: return "triangle-midpass";
			case SelectivityFilter::TriangleHighpass: return "triangle-highpass";
			case SelectivityFilter::HomographLowpass: return "homograph-lowpass";
		} 
	return "invalid";
	}

	bool operator==(SelectivityFilter::Enumeration ref) { return value == ref; }
	SelectivityFilter& operator=(SelectivityFilter::Enumeration ref) { value = ref; return *this; }

protected:
	Enumeration value;
};


class SectorCopy
{
public:
	enum Enumeration {Full, Row, None};
	
	SectorCopy() {} 
	SectorCopy(std::string s) {
		if (s == "full") { value = SectorCopy::Full; }
		else if (s == "row") { value = SectorCopy::Row; }
		else if (s == "none") { value = SectorCopy::None; }
	}

	std::string toString() {
		switch (value) {
			case SectorCopy::Full: return "full";
			case SectorCopy::Row: return "row";
			case SectorCopy::None: return "none";
		}
		return "invalid";
	}

	bool operator==(SectorCopy::Enumeration ref) { return value == ref; }
	SectorCopy& operator=(SectorCopy::Enumeration ref) { value = ref; return *this; }

protected:
	enum Enumeration value;
};



class RefinementType 
{
public:
	enum Enumeration {Random, Column, Tree};

	RefinementType() {}
	RefinementType(std::string s) {
		if (s == "random") { value = Random; }
		else if (s == "column") { value = Column; }
		else if (s == "tree") { value = Tree; }
	}

	std::string toString() {
		switch (value) {
			case Random: return "random";
			case Column: return "column";
			case Tree: return "tree";
		}
		return "invalid";
	}

	bool operator==(RefinementType::Enumeration ref) { return value == ref; }
	RefinementType& operator=(RefinementType::Enumeration ref) { value = ref; return *this; }


protected:
	Enumeration value;
};