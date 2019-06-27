#ifndef _MSA_H
#define _MSA_H
#include "MSADef.h"
#include "MSAGuideTree.h"

#include "SafeVector.h"
#include "MultiSequence.h"
#include "ScoreType.h"
#include "ProbabilisticModel.h"
#include "SparseMatrix.h"
#include <string>
#include "AlignGraph.h"

using namespace std;

class MSAGuideTree;
struct TreeNode;
class MSA{
public:
	MSA(int argc, char* argv[]);
	~MSA();

	static void getSysTime(double * dtime);
	MSAGuideTree* getGuideTree() {
		return tree;
	}
	int * getSeqsWeights() {
		return seqsWeights;
	}
private:
	//print usage
	void printUsage();
	//do multiple sequence alignment
	void doAlign();

	//for sequence weights
	void createSeqsWeights(int seqsNum);
	void releaseSeqsWeights();

	//weights of sequences
	int * seqsWeights;
	//guide tree for progressive alignment
	MSAGuideTree* tree;
	//alignment graph for non-progressive alignment
	AlignGraph *graph;
	//output file
	string alignOutFileName;
	std::ostream* alignOutFile;
private:
	SafeVector<string> ParseParams(int argc, char *argv[]);
	void PrintParameters(const char *message, const VF &initDistrib,
			const VF &gapOpen, const VF &gapExtend, const VVF &emitPairs,
			const VF &emitSingle, const char *filename);

	SafeVector<string> PostProbsParseParams(int argc, char **argv);

	void ReadParameters();

    int ModelAdjustmentTest( MultiSequence *sequences );//Determine the Model
		string Alter_ModelAdjustmentTest( MultiSequence *sequences );//Determine the Model

    MultiSequence *pdoAlign(MultiSequence *sequencen, const int variance_mean); //progressive alignment
    MultiSequence *npdoAlign(MultiSequence *sequencen, const int variance_mean); //non-progressive alignment


    //consistency transformation
	SafeVector<SafeVector<SparseMatrix *> > DoRelaxation( MultiSequence *sequences,
			SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices);

	void Relax(SparseMatrix *matXZ, SparseMatrix *matZY,VF &posterior);//unweight
	void Relax1(SparseMatrix *matXZ, SparseMatrix *matZY,VF &posterior);//unweight
	MultiSequence *AlignAlignments(MultiSequence *align1, MultiSequence *align2,
			const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
			const ProbabilisticModel &model, float &alignscore, bool nflag);

	//progreesive functions
	MultiSequence* ProcessTree(TreeNode *tree, MultiSequence *sequences,
			const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
			const ProbabilisticModel &model);
	MultiSequence* ComputeFinalAlignment(MSAGuideTree*tree, MultiSequence *sequences,
		const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
		const ProbabilisticModel &model, const int pid);
	int DoIterativeRefinement(
			const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
			const ProbabilisticModel &model, MultiSequence* &alignment);

	//non-progressive
	void ArrangePosteriorProbs(MultiSequence *sequences,
		const ProbabilisticModel &model,
		SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
		VVF &distances, const int variance_mean);//Compute HMMs
	MultiSequence* ComputeGraph(MultiSequence *sequences, const ProbabilisticModel &model,
		SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, VVF &distances); //Compute Alignment Graph
	void DoRefinement(MultiSequence* &alignment, SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
		 const ProbabilisticModel &model, VVF distances);
	void FindSimilar(VVF distances, vector<set<int> > &SimSeqs);


	void WriteAnnotation(MultiSequence *alignment,
			const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices);
	int ComputeScore(const SafeVector<pair<int, int> > &active,
			const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices);

	bool GetInteger(char *data, int *val);
	bool GetFloat(char *data, float *val);

#ifdef _OPENMP
	//private struct
	struct SeqsPair {
		int seq1;
		int seq2;
	};
	int numPairs;
	SeqsPair* seqsPairs;
#endif
};

#endif
