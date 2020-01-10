#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <list>
#include <set>
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <iomanip>
#include "MSA.h"
#include "MSAClusterTree.h"
#include "Defaults.h"
#define MAX_ARR 10000
#include <sys/time.h>
#include <time.h>
using namespace std;
#ifdef _OPENMP
#include <omp.h>
#endif

string parametersInputFilename = "";
string parametersOutputFilename = "no training";
string annotationFilename = "";
bool use_pnpprobs_big = true;
bool enableVerbose = false;
bool enableRunningtime = true;
bool enableAnnotation = false;
bool enableClustalWOutput = false;
bool enableAlignOrder = false;
int numConsistencyReps = 2;
int numPreTrainingReps = 0;
int numIterativeRefinementReps = 100;
bool just_getpid = false;
float cutoff = 0;

///////////////////////////////
// local pair-HMM and double affine pair-HMM variables
//////////////////////////////
VF initDistrib(NumMatrixTypes);
VF gapOpen(2 * NumInsertStates);
VF gapExtend(2 * NumInsertStates);
VVF emitPairs(256, VF(256, 1e-10));//subsititution matrix for double affine pair-HMM
VF emitSingle(256, 1e-5);//subsititution matrix for double affine pair-HMM

string alphabet = alphabetDefault;

const int MIN_PRETRAINING_REPS = 0;
const int MAX_PRETRAINING_REPS = 20;
const int MIN_CONSISTENCY_REPS = 0;
const int MAX_CONSISTENCY_REPS = 5;
const int MIN_ITERATIVE_REFINEMENT_REPS = 0;
const int MAX_ITERATIVE_REFINEMENT_REPS = 1000;

string posteriorProbsFilename = "";
bool allscores = true;
string infilename;

///////////////////////////////
// global pair-HMM variables
//////////////////////////////
int flag_gui = 0;   //0: no gui related o/p
//1: gui related o/p generated
int flag_ppscore = 0; //0: no pp score sequence added to o/p fasta alignment
//1: pp score seq added to o/p fasta alignment

float g_gap_open1, g_gap_open2, g_gap_ext1, g_gap_ext2;
char *aminos, *bases, matrixtype[20] = "gonnet_160";
int subst_index[26];

double sub_matrix[26][26];      //subsititution matrix for global pair-HMM
double normalized_matrix[26][26];// be used to compute sequences' similarity score
int firstread = 0;		//this makes sure that matrices are read only once

float TEMPERATURE = 5;
int MATRIXTYPE = 160;
int prot_nuc = 0;		//0=prot, 1=nucleotide

float GAPOPEN = 0;
float GAPEXT = 0;
int numThreads = 0;

double startTime = 0;
double lastTime = 0;
double timeUsed = 0;

//argument support
typedef struct {
	char input[30];
	int matrix;
	int N;
	float T;
	float beta;
	char opt;			//can be 'P' or 'M'
	float gapopen;
	float gapext;
} argument_decl;

argument_decl argument;

extern inline void read_sustitution_matrix(char *fileName);
extern void setmatrixtype(int le);
extern inline int matrixtype_to_int();
extern inline void read_dna_matrix();
extern inline void read_vtml_la_matrix();
extern void init_arguments();

unsigned long int GetTime()
{
    struct timeval tv;
    gettimeofday ( &tv, NULL );
    return tv.tv_usec + ( tv.tv_sec * 1000000L );
};

double GetElapsedTime ( unsigned long int startTime )
{
    return ( GetTime() - startTime ) / 1000000.0L;
};

MSA::MSA(int argc, char* argv[]) {
	//parse program parameters
	SafeVector<string> sequenceNames = ParseParams(argc, argv);

	//initialize arguments for partition function
	init_arguments();

	ReadParameters();

	//read the input sequences
	MultiSequence *sequences = new MultiSequence();
	assert(sequences);
	for (int i = 0; i < (int) sequenceNames.size(); i++) {
		sequences->LoadMFA(sequenceNames[i], true);
	}
	//allocate space for sequence weights
	this->seqsWeights = new int[sequences->GetNumSequences()];
	//initilaize the guide tree
	this->tree = NULL;
	//initilaize the alignment graph
	this->graph = NULL;

	//initilaize parameters for OPENMP
#ifdef _OPENMP
	if(numThreads <= 0) {
		numThreads = omp_get_num_procs();
	}
	//set OpenMP to use dynamic number of threads which is equal to the number of processor cores on the host
	omp_set_num_threads(numThreads);
#endif
//getPID !!!!
if (just_getpid){
	string pid = Alter_ModelAdjustmentTest( sequences, 1.0 );
	/*
	ofstream outfile;
	outfile.open("./tmp/tmp_pid.txt", ios::out);
	// 向文件写入用户输入的数据
	outfile << pid << endl;
	// 关闭打开的文件
	outfile.close();
*/
	cout << pid << endl;
	return ;
}


	//model determination
	startTime = GetTime();
	int variance_mean = ModelAdjustmentTest( sequences );
	timeUsed = GetElapsedTime ( startTime );
	// now, we can perform the alignments and write them out
	MultiSequence *alignment; //variance_mean%10 > 0

  if(use_pnpprobs_big){
    alignment = pdoAlign( sequences, variance_mean ); //progressive
  }
	else alignment = npdoAlign( sequences, variance_mean ); //non-progressive

	//write the alignment results to standard output
	alignment->WriteMFA(*alignOutFile);
	//release resources
	delete alignment;
	delete[] this->seqsWeights;
	delete sequences;
}

MSA::~MSA() {
	/*close the output file*/
	if (alignOutFileName.length() > 0) {
		((std::ofstream*) alignOutFile)->close();
	}
}

/////////////////////////////////////////////////////////////////
// ParseParams()
//
// Parse all command-line options.
/////////////////////////////////////////////////////////////////

void MSA::printUsage() {
	cerr
			<< "************************************************************************"
			<< endl
			<< "This program is modified from PNPProbs by M.M. Kuang (mmkuang@cs.hku.hk)."
			<< endl
			<< "PNPProbs is a multiple sequence alignment program"
			<< endl
			<< "combining progressive and non-progressive alignment methods."
			<< endl
			<< "If any comments or problems, please contact: Yongtao YE (ytye@cs.hku.hk)."
			<< endl
			<< "*************************************************************************"
			<< endl << "Usage:" << endl
			<< "       pnpprobs [infile] ...[OPTION]... " << endl << endl
			<< "Description:" << endl
			<< "       Align sequences in multi-FASTA format" << endl << endl
			<< "       -o, --outfile <string>" << endl
			<< "              specify the output file name (STDOUT by default)" << endl
			<< "       -p, --program <string>" << endl
			<< "              choose progressive or no-progressive alignment parts (progressive by default)"
			<< endl << "       -clustalw" << endl
			<< "              use CLUSTALW output format instead of FASTA format"
			<< endl << endl << "       -c, --consistency REPS" << endl
			<< "              use " << MIN_CONSISTENCY_REPS << " <= REPS <= "
			<< MAX_CONSISTENCY_REPS << " (default: " << numConsistencyReps
			<< ") passes of consistency transformation" << endl << endl
			<< "       -ir, --iterative-refinement REPS" << endl
			<< "              use " << MIN_ITERATIVE_REFINEMENT_REPS
			<< " <= REPS <= " << MAX_ITERATIVE_REFINEMENT_REPS << " (default: "
			<< numIterativeRefinementReps << ") passes of iterative-refinement"

			<< endl << endl << "       -timeon, -timeoff" << endl
			<< "              report program running times of each step (default:on) "
			<< endl << endl << "       -v, --verbose" << endl
			<< "              report progress while aligning (default: "
			<< (enableVerbose ? "on" : "off") << ")" << endl << endl
			<< "       -annot FILENAME" << endl
			<< "              write annotation for multiple alignment to FILENAME"
			<< endl << endl << "       -a, --alignment-order" << endl
			<< "              print sequences in alignment order rather than input order (default: "
			<< (enableAlignOrder ? "on" : "off") << ")" << endl
			<< "       -version " << endl
			<< "              print out version of PNPProbs " << endl << endl;
}

SafeVector<string> MSA::ParseParams(int argc, char **argv) {
	if (argc < 2) {
		printUsage();
		exit(1);
	}
	SafeVector<string> sequenceNames;
	int tempInt;
	float tempFloat;

	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			//help
			if (!strcmp(argv[i], "-help") || !strcmp(argv[i], "-?")) {
				printUsage();
				exit(1);
			//output file name
			} else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--outfile")) {
				if (i < argc - 1) {
					alignOutFileName = argv[++i];	//get the file name
				} else {
					cerr << "ERROR: String expected for option " << argv[i]
							<< endl;
					exit(1);
				}
        //使用哪个程序
        }else if (!strcmp (argv[i], "-p") || !strcmp (argv[i], "--program")){
                  if (!GetInteger(argv[++i], &tempInt)) {
              cerr << "ERROR: Invalid integer following option "
                  << argv[i - 1] << ": " << argv[i] << endl;
              exit(1);
          } else {
            if (tempInt > 1 || tempInt < 0) {
              cerr << "ERROR: For option " << argv[i - 1]
                << ", integer must be 0 or 1." << endl;
                exit(1);
            } else {
              if(tempInt == 0){
                use_pnpprobs_big = true;
              }else{
                use_pnpprobs_big = false;
              }
            }
          }

					//getpid
				}else if (!strcmp (argv[i], "-G") || !strcmp (argv[i], "--getPID")){
										just_getpid = true;
			// number of consistency transformations
			} else if (!strcmp(argv[i], "-c")
					|| !strcmp(argv[i], "--consistency")) {
				if (i < argc - 1) {
					if (!GetInteger(argv[++i], &tempInt)) {
						cerr << "ERROR: Invalid integer following option "
								<< argv[i - 1] << ": " << argv[i] << endl;
						exit(1);
					} else {
						if (tempInt < MIN_CONSISTENCY_REPS
								|| tempInt > MAX_CONSISTENCY_REPS) {
							cerr << "ERROR: For option " << argv[i - 1]
									<< ", integer must be between "
									<< MIN_CONSISTENCY_REPS << " and "
									<< MAX_CONSISTENCY_REPS << "." << endl;
							exit(1);
						} else {
							numConsistencyReps = tempInt;
						}
					}
				} else {
					cerr << "ERROR: Integer expected for option " << argv[i]
							<< endl;
					exit(1);
				}
			}

			// number of randomized partitioning iterative refinement passes
			else if (!strcmp(argv[i], "-ir")
					|| !strcmp(argv[i], "--iterative-refinement")) {
				if (i < argc - 1) {
					if (!GetInteger(argv[++i], &tempInt)) {
						cerr << "ERROR: Invalid integer following option "
								<< argv[i - 1] << ": " << argv[i] << endl;
						exit(1);
					} else {
						if (tempInt < MIN_ITERATIVE_REFINEMENT_REPS
								|| tempInt > MAX_ITERATIVE_REFINEMENT_REPS) {
							cerr << "ERROR: For option " << argv[i - 1]
									<< ", integer must be between "
									<< MIN_ITERATIVE_REFINEMENT_REPS << " and "
									<< MAX_ITERATIVE_REFINEMENT_REPS << "."
									<< endl;
							exit(1);
						} else
							numIterativeRefinementReps = tempInt;
					}
				} else {
					cerr << "ERROR: Integer expected for option " << argv[i]
							<< endl;
					exit(1);
				}
			}

			// annotation files
			else if (!strcmp(argv[i], "-annot")) {
				enableAnnotation = true;
				if (i < argc - 1) {
					annotationFilename = argv[++i];
				} else {
					cerr << "ERROR: FILENAME expected for option " << argv[i]
							<< endl;
					exit(1);
				}
			}

			// clustalw output format
			else if (!strcmp(argv[i], "-clustalw")) {
				enableClustalWOutput = true;
			}

			// cutoff
			else if (!strcmp(argv[i], "-co") || !strcmp(argv[i], "--cutoff")) {
				if (i < argc - 1) {
					if (!GetFloat(argv[++i], &tempFloat)) {
						cerr << "ERROR: Invalid floating-point value following option "
								<< argv[i - 1] << ": " << argv[i] << endl;
						exit(1);
					} else {
						if (tempFloat < 0 || tempFloat > 1) {
							cerr << "ERROR: For option " << argv[i - 1]
									<< ", floating-point value must be between 0 and 1."
									<< endl;
							exit(1);
						} else
							cutoff = tempFloat;
					}
				} else {
					cerr << "ERROR: Floating-point value expected for option "
							<< argv[i] << endl;
					exit(1);
				}
			}

			// verbose reporting
			else if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) {
				enableVerbose = true;
			}

			// runnint time reporting
			else if (!strcmp(argv[i], "-timeoff")) {
				enableRunningtime = false;
			}
			else if (!strcmp(argv[i], "-timeon")) {
				enableRunningtime = true;
			}

			// alignment order
			else if (!strcmp(argv[i], "-a")
					|| !strcmp(argv[i], "--alignment-order")) {
				enableAlignOrder = true;
			}

			//print out version
			else if (!strcmp(argv[i], "-version")) {
				cerr << "PNPProbs version " << VERSION << endl;
				exit(1);
			}
			// bad arguments
			else {
				cerr << "ERROR: Unrecognized option: " << argv[i] << endl;
				exit(1);
			}
		} else {
			sequenceNames.push_back(string(argv[i]));
		}
	}

	/*check the output file name*/
	//cerr << "-------------------------------------" << endl;
	if (alignOutFileName.length() == 0) {
		//cerr << "The final alignments will be printed out to STDOUT" << endl;
		alignOutFile = &std::cout;
	} else {
		//cerr << "Open the output file " << alignOutFileName << endl;
		alignOutFile = new ofstream(alignOutFileName.c_str(),
				ios::binary | ios::out | ios::trunc);
	}
	//cerr << "-------------------------------------" << endl;
	return sequenceNames;
}

/////////////////////////////////////////////////////////////////
// ReadParameters()
//
// Read initial distribution, transition, and emission
// parameters from a file.
/////////////////////////////////////////////////////////////////

void MSA::ReadParameters() {

	ifstream data;

	emitPairs = VVF(256, VF(256, 1e-10));
	emitSingle = VF(256, 1e-5);

	// read initial state distribution and transition parameters
	if (parametersInputFilename == string("")) {
		if (NumInsertStates == 1) {
			for (int i = 0; i < NumMatrixTypes; i++)
				initDistrib[i] = initDistrib1Default[i];
			for (int i = 0; i < 2 * NumInsertStates; i++)
				gapOpen[i] = gapOpen1Default[i];
			for (int i = 0; i < 2 * NumInsertStates; i++)
				gapExtend[i] = gapExtend1Default[i];
		} else if (NumInsertStates == 2) {
			for (int i = 0; i < NumMatrixTypes; i++)
				initDistrib[i] = initDistrib2Default[i];
			for (int i = 0; i < 2 * NumInsertStates; i++)
				gapOpen[i] = gapOpen2Default[i];
			for (int i = 0; i < 2 * NumInsertStates; i++)
				gapExtend[i] = gapExtend2Default[i];
		} else {
			cerr
					<< "ERROR: No default initial distribution/parameter settings exist"
					<< endl << "       for " << NumInsertStates
					<< " pairs of insert states.  Use --paramfile." << endl;
			exit(1);
		}

		alphabet = alphabetDefault;

		for (int i = 0; i < (int) alphabet.length(); i++) {
			emitSingle[(unsigned char) tolower(alphabet[i])] =
					emitSingleDefault[i];
			emitSingle[(unsigned char) toupper(alphabet[i])] =
					emitSingleDefault[i];
			for (int j = 0; j <= i; j++) {
				emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) tolower(
						alphabet[j])] = emitPairsDefault[i][j];
				emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) toupper(
						alphabet[j])] = emitPairsDefault[i][j];
				emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) tolower(
						alphabet[j])] = emitPairsDefault[i][j];
				emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) toupper(
						alphabet[j])] = emitPairsDefault[i][j];
				emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) tolower(
						alphabet[i])] = emitPairsDefault[i][j];
				emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) toupper(
						alphabet[i])] = emitPairsDefault[i][j];
				emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) tolower(
						alphabet[i])] = emitPairsDefault[i][j];
				emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) toupper(
						alphabet[i])] = emitPairsDefault[i][j];
			}
		}
	} else {
		data.open(parametersInputFilename.c_str());
		if (data.fail()) {
			cerr << "ERROR: Unable to read parameter file: "
					<< parametersInputFilename << endl;
			exit(1);
		}

		string line[3];
		for (int i = 0; i < 3; i++) {
			if (!getline(data, line[i])) {
				cerr
						<< "ERROR: Unable to read transition parameters from parameter file: "
						<< parametersInputFilename << endl;
				exit(1);
			}
		}
		istringstream data2;
		data2.clear();
		data2.str(line[0]);
		for (int i = 0; i < NumMatrixTypes; i++)
			data2 >> initDistrib[i];
		data2.clear();
		data2.str(line[1]);
		for (int i = 0; i < 2 * NumInsertStates; i++)
			data2 >> gapOpen[i];
		data2.clear();
		data2.str(line[2]);
		for (int i = 0; i < 2 * NumInsertStates; i++)
			data2 >> gapExtend[i];

		if (!getline(data, line[0])) {
			cerr << "ERROR: Unable to read alphabet from scoring matrix file: "
					<< parametersInputFilename << endl;
			exit(1);
		}

		// read alphabet as concatenation of all characters on alphabet line
		alphabet = "";
		string token;
		data2.clear();
		data2.str(line[0]);
		while (data2 >> token)
			alphabet += token;

		for (int i = 0; i < (int) alphabet.size(); i++) {
			for (int j = 0; j <= i; j++) {
				float val;
				data >> val;
				emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) tolower(
						alphabet[j])] = val;
				emitPairs[(unsigned char) tolower(alphabet[i])][(unsigned char) toupper(
						alphabet[j])] = val;
				emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) tolower(
						alphabet[j])] = val;
				emitPairs[(unsigned char) toupper(alphabet[i])][(unsigned char) toupper(
						alphabet[j])] = val;
				emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) tolower(
						alphabet[i])] = val;
				emitPairs[(unsigned char) tolower(alphabet[j])][(unsigned char) toupper(
						alphabet[i])] = val;
				emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) tolower(
						alphabet[i])] = val;
				emitPairs[(unsigned char) toupper(alphabet[j])][(unsigned char) toupper(
						alphabet[i])] = val;
			}
		}

		for (int i = 0; i < (int) alphabet.size(); i++) {
			float val;
			data >> val;
			emitSingle[(unsigned char) tolower(alphabet[i])] = val;
			emitSingle[(unsigned char) toupper(alphabet[i])] = val;
		}
		data.close();
	}
}



/////////////////////////////////////////////////////////////////
// PrintParameters()
//
// Prints MSAPROBS parameters to STDERR.  If a filename is
// specified, then the parameters are also written to the file.
/////////////////////////////////////////////////////////////////

void MSA::PrintParameters(const char *message, const VF &initDistrib,
		const VF &gapOpen, const VF &gapExtend, const VVF &emitPairs,
		const VF &emitSingle, const char *filename) {

	// print parameters to the screen
	cerr << message << endl << "    initDistrib[] = { ";
	for (int i = 0; i < NumMatrixTypes; i++)
		cerr << setprecision(10) << initDistrib[i] << " ";
	cerr << "}" << endl << "        gapOpen[] = { ";
	for (int i = 0; i < NumInsertStates * 2; i++)
		cerr << setprecision(10) << gapOpen[i] << " ";
	cerr << "}" << endl << "      gapExtend[] = { ";
	for (int i = 0; i < NumInsertStates * 2; i++)
		cerr << setprecision(10) << gapExtend[i] << " ";
	cerr << "}" << endl << endl;

	// if a file name is specified
	if (filename) {

		// attempt to open the file for writing
		FILE *file = fopen(filename, "w");
		if (!file) {
			cerr << "ERROR: Unable to write parameter file: " << filename
					<< endl;
			exit(1);
		}

		// if successful, then write the parameters to the file
		for (int i = 0; i < NumMatrixTypes; i++)
			fprintf(file, "%.10f ", initDistrib[i]);
		fprintf(file, "\n");
		for (int i = 0; i < 2 * NumInsertStates; i++)
			fprintf(file, "%.10f ", gapOpen[i]);
		fprintf(file, "\n");
		for (int i = 0; i < 2 * NumInsertStates; i++)
			fprintf(file, "%.10f ", gapExtend[i]);
		fprintf(file, "\n");
		fprintf(file, "%s\n", alphabet.c_str());
		for (int i = 0; i < (int) alphabet.size(); i++) {
			for (int j = 0; j <= i; j++)
				fprintf(file, "%.10f ",
						emitPairs[(unsigned char) alphabet[i]][(unsigned char) alphabet[j]]);
			fprintf(file, "\n");
		}
		for (int i = 0; i < (int) alphabet.size(); i++)
			fprintf(file, "%.10f ", emitSingle[(unsigned char) alphabet[i]]);
		fprintf(file, "\n");
		fclose(file);
	}
}


/////////////////////////////////////////////////////////////////
// Alter_ModelAdjustmentTest ()
//
// Computes the average percent identity for a particular family.
/////////////////////////////////////////////////////////////////

string MSA::Alter_ModelAdjustmentTest(MultiSequence *sequences, float theta){
	assert(sequences);

	const int numSeqs = sequences->GetNumSequences();
	ProbabilisticModel model(initDistrib, gapOpen, gapExtend, emitPairs,emitSingle);
    float identity = 0;
		int max_length_pair = 0;
		int tmp_sp_idx = 0;
		float tmp_sp = 0;
		float* new_final_arr = new float[MAX_ARR]();
#ifdef _OPENMP
	int pairIdx = 0;
	numPairs = (numSeqs - 1) * numSeqs / 2;
	seqsPairs = new SeqsPair[numPairs];
	for(int a = 0; a < numSeqs; a++) {
		for(int b = a + 1; b < numSeqs; b++) {
			seqsPairs[pairIdx].seq1 = a;
			seqsPairs[pairIdx].seq2 = b;
			pairIdx++;
		}
	}
#endif

	float* PIDs = new float[ (numSeqs - 1) * numSeqs / 2 * sizeof(float) ];
 float* SOPs = new float[ (numSeqs - 1) * numSeqs / 2 * sizeof(float) ];
 int avg_length = 0;
#ifdef _OPENMP
#pragma omp parallel for private(pairIdx) default(shared) schedule(dynamic)
	for(pairIdx = 0; pairIdx < numPairs; pairIdx++) {
		int a= seqsPairs[pairIdx].seq1;
		int b = seqsPairs[pairIdx].seq2;
		if(enableVerbose) {
#pragma omp critical
			cerr <<"tid "<<omp_get_thread_num()<<" a "<<a<<" b "<<b<<endl;
		}
#else
	int pairIdx = -1;
	for (int a = 0; a < numSeqs - 1; a++) {
		for (int b = a + 1; b < numSeqs; b++) {
			pairIdx++;
#endif
 			int num_idx = 0;
			Sequence *seq1 = sequences->GetSequence(a);
			Sequence *seq2 = sequences->GetSequence(b);
			pair<SafeVector<char> *, float> alignment = model.ComputeViterbiAlignment(seq1,seq2);
			SafeVector<char>::iterator iter1 = seq1->GetDataPtr();
    		SafeVector<char>::iterator iter2 = seq2->GetDataPtr();
				avg_length += alignment.first -> size();
				if(alignment.first -> size() > max_length_pair){
					max_length_pair = alignment.first -> size();
				}
            float N_correct_match = 0;
						float N_emit = 0;
            int i = 1;int j = 1;
			for (SafeVector<char>::iterator iter = alignment.first->begin();
				iter != alignment.first->end(); ++iter){
				if (*iter == 'B'){
					unsigned char c1 = (unsigned char) iter1[i++];
					unsigned char c2 = (unsigned char) iter2[j++];
   					if(c1==c2) N_correct_match += 1;
						N_emit += emitPairsDefault[alphabetDefault.find(c1)][alphabetDefault.find(c2)];
						if (BLOSUM62[alphabetDefault.find(c1)][alphabetDefault.find(c2)] < 10){
							new_final_arr[num_idx] += BLOSUM62[alphabetDefault.find(c1)][alphabetDefault.find(c2)];
							tmp_sp += BLOSUM62[alphabetDefault.find(c1)][alphabetDefault.find(c2)];
						}else{
							new_final_arr[num_idx] += 0;
						}

						tmp_sp_idx += 1;
						num_idx ++;
				}

				else if(*iter == 'X'){
          i++;
          num_idx ++;
					tmp_sp_idx += 1;
        }
 				else if(*iter == 'Y'){
           j++;
           num_idx ++;
					 tmp_sp_idx += 1;
        }
            }
            if(i!= seq1->GetLength()+1 || j!= seq2->GetLength() + 1 ) cerr << "percent identity error"<< endl;
            PIDs[pairIdx] = N_correct_match / alignment.first->size();
            identity += N_correct_match / alignment.first->size();
						SOPs[pairIdx] =  N_emit;
			delete alignment.first;
#ifndef _OPENMP
		}
#endif
	}

tmp_sp /= tmp_sp_idx;
	identity /=  ( (numSeqs-1)*numSeqs/2 );//average percent identity
	avg_length /= ( (numSeqs-1)*numSeqs/2 );
float peak_length_ = 0;
 int tmp_i = 0;
 for(tmp_i;  tmp_i < max_length_pair; tmp_i ++){
	 new_final_arr[tmp_i] /= ((numSeqs-1)*numSeqs/2 );
	 if (theta <= new_final_arr[tmp_i]){
		 peak_length_ += 1;
	 }
 }
 peak_length_ /= max_length_pair;

	float variance = 0;

	for (int k = 0; k < (numSeqs-1)*numSeqs/2; k++) {
		variance += (PIDs[k]-identity)*(PIDs[k]-identity);
	}
	variance /= ( (numSeqs-1)*numSeqs/2 );
	variance = sqrt(variance);
    float factor = 2 * (float)numSeqs - (float) avg_length;

    return to_string(identity) + "\t" + to_string(variance) + "\t" + to_string(numSeqs) + "\t" + to_string(avg_length) + "\t" + to_string(tmp_sp) + "\t" + to_string(peak_length_) + "\t" + to_string(factor);
}


/////////////////////////////////////////////////////////////////
// ModelAdjustmentTest ()
//
// Computes the average percent identity for a particular family.
// Divergent  (<=25%) return 0(0-20%), 1(20-25%)
// Medium   (25%-40%) return 2
// Similar  (40%-70%) return 3
// High Similar(>70%) return 4
/////////////////////////////////////////////////////////////////

int MSA::ModelAdjustmentTest(MultiSequence *sequences){
	assert(sequences);

	//get the number of sequences
	const int numSeqs = sequences->GetNumSequences();
	//initialize hmm model
	ProbabilisticModel model(initDistrib, gapOpen, gapExtend, emitPairs,emitSingle);
    //average identity for all sequences
    float identity = 0;

#ifdef _OPENMP
	//calculate sequence pairs for openmp model
	int pairIdx = 0;
	numPairs = (numSeqs - 1) * numSeqs / 2;
	seqsPairs = new SeqsPair[numPairs];
	for(int a = 0; a < numSeqs; a++) {
		for(int b = a + 1; b < numSeqs; b++) {
			seqsPairs[pairIdx].seq1 = a;
			seqsPairs[pairIdx].seq2 = b;
			pairIdx++;
		}
	}
#endif

	// do all pairwise alignments for family similarity
	//store percent identity of every pair sequences
	float* PIDs = new float[ (numSeqs - 1) * numSeqs / 2 * sizeof(float) ];

#ifdef _OPENMP
#pragma omp parallel for private(pairIdx) default(shared) schedule(dynamic)
	for(pairIdx = 0; pairIdx < numPairs; pairIdx++) {
		int a= seqsPairs[pairIdx].seq1;
		int b = seqsPairs[pairIdx].seq2;
		if(enableVerbose) {
#pragma omp critical
			cerr <<"tid "<<omp_get_thread_num()<<" a "<<a<<" b "<<b<<endl;
		}
#else
	int pairIdx = -1;
	for (int a = 0; a < numSeqs - 1; a++) {
		for (int b = a + 1; b < numSeqs; b++) {
			pairIdx++;
#endif
			Sequence *seq1 = sequences->GetSequence(a);
			Sequence *seq2 = sequences->GetSequence(b);
			pair<SafeVector<char> *, float> alignment = model.ComputeViterbiAlignment(seq1,seq2);
			SafeVector<char>::iterator iter1 = seq1->GetDataPtr();
    		SafeVector<char>::iterator iter2 = seq2->GetDataPtr();

            float N_correct_match = 0;
			//float N_alignment = 0;
            int i = 1;int j = 1;
			for (SafeVector<char>::iterator iter = alignment.first->begin();
				iter != alignment.first->end(); ++iter){
				//N_alignment += 1;
				if (*iter == 'B'){
					unsigned char c1 = (unsigned char) iter1[i++];
					unsigned char c2 = (unsigned char) iter2[j++];
   					if(c1==c2) N_correct_match += 1;
				}
                else if(*iter == 'X') i++;
 				else if(*iter == 'Y') j++;
            }
            if(i!= seq1->GetLength()+1 || j!= seq2->GetLength() + 1 ) cerr << "percent identity error"<< endl;
            PIDs[pairIdx] = N_correct_match / alignment.first->size();
            identity += N_correct_match / alignment.first->size();
			//identity += N_correct_match / alignment.first->size();
			delete alignment.first;
#ifndef _OPENMP
		}
#endif
	}
	identity /=  ( (numSeqs-1)*numSeqs/2 );//average percent identity

	//compute the variance
	float variance = 0;

	for (int k = 0; k < (numSeqs-1)*numSeqs/2; k++) {
		variance += (PIDs[k]-identity)*(PIDs[k]-identity);
	}
	variance /= ( (numSeqs-1)*numSeqs/2 );
	variance = sqrt(variance);

    //adjust the parameter of leaving RX/RY in random pair-HMM


    if( identity <= 0.125 ) initDistrib[2] = 0.108854;//0.173854;
    else if( identity <= 0.15 ) initDistrib[2] = 0.132548;
    else if( identity <= 0.175 ) initDistrib[2] = 0.165248;//0.167248;//0.161748;//0.191948;
	else if( identity <= 0.2 ) initDistrib[2] = 0.168284;
	else if( identity <= 0.25 ) initDistrib[2] = 0.170705;
	else if( identity <= 0.3 ) initDistrib[2] = 0.100675;
	else if( identity <= 0.35 ) initDistrib[2] = 0.090755;
	else if( identity <= 0.4 ) initDistrib[2] = 0.146188;
    else if( identity <= 0.45 ) initDistrib[2] = 0.167858;
	else if( identity <= 0.5) initDistrib[2] = 0.250769;


    int variance_mean = 0;
//	if( variance > 0.028 ) variance_mean = 10;
    if( variance > 0.115 ) variance_mean = 10;

    if( identity <= 0.18 ) return variance_mean + 0;
    else if( identity <= 0.25 ) return variance_mean + 1;
    else if( identity <= 0.4) return variance_mean + 2;
	else if( identity <= 0.7) return variance_mean + 3;
    else return variance_mean + 4;
}

/////////////////////////////////////////////////////////////////
// doAlign()
//
// First computes all pairwise posterior probability matrices.
// Then, computes new parameters if training, or a final
// alignment, otherwise.
/////////////////////////////////////////////////////////////////

extern VF *ComputePostProbs(int a, int b, string seq1, string seq2);

//progressive alignment
MultiSequence* MSA::pdoAlign( MultiSequence *sequences, const int variance_mean ) {
	assert(sequences);

	//get the number of sequences
	const int numSeqs = sequences->GetNumSequences();
	//creat sparseMatrices
    SafeVector<SafeVector<SparseMatrix *> > sparseMatrices(numSeqs,
			SafeVector<SparseMatrix *>(numSeqs, NULL));
    //create distance matrix
 	VVF distances(numSeqs, VF(numSeqs, 0));

 	startTime = GetTime();
#ifdef _OPENMP
	//calculate sequence pairs for openmp model
	int pairIdx = 0;
	numPairs = (numSeqs - 1) * numSeqs / 2;
	seqsPairs = new SeqsPair[numPairs];
	for(int a = 0; a < numSeqs; a++) {
		for(int b = a + 1; b < numSeqs; b++) {
			seqsPairs[pairIdx].seq1 = a;
			seqsPairs[pairIdx].seq2 = b;
			pairIdx++;
		}
	}
#endif

 	ProbabilisticModel model(initDistrib, gapOpen, gapExtend, emitPairs,emitSingle);
 	int pid = variance_mean%10;
 	int vpid = variance_mean/10;

	// do all pairwise alignments for posterior probability matrices
#ifdef _OPENMP
#pragma omp parallel for private(pairIdx) default(shared) schedule(dynamic)
	for(pairIdx = 0; pairIdx < numPairs; pairIdx++) {
		int a= seqsPairs[pairIdx].seq1;
		int b = seqsPairs[pairIdx].seq2;
		if(enableVerbose) {
#pragma omp critical
			cerr <<"tid "<<omp_get_thread_num()<<" a "<<a<<" b "<<b<<endl;
		}
#else
	for (int a = 0; a < numSeqs - 1; a++) {
		for (int b = a + 1; b < numSeqs; b++) {
#endif
			Sequence *seq1 = sequences->GetSequence(a);
			Sequence *seq2 = sequences->GetSequence(b);

			//posterior probability matrix
			VF* posterior;

//low similarity use local pair-HMM
			if(pid == 2){
				// compute forward and backward probabilities
				VF *forward = model.ComputeForwardMatrix(seq1, seq2,false);
				assert(forward);
				VF *backward = model.ComputeBackwardMatrix(seq1, seq2,false);
				assert(backward);
				// compute posterior probability
				posterior = model.ComputePosteriorMatrix(seq1, seq2, *forward,*backward, false);
   				delete forward;
				delete backward;

			}
//high similarity use global pair-HMM
			else if(pid >= 3) posterior = ::ComputePostProbs(a, b, seq1->GetString(),seq2->GetString());

//extreme low use combined model
			else{

//double affine pair-HMM
				// compute forward and backward probabilities
				VF *forward = model.ComputeForwardMatrix(seq1, seq2);
				assert(forward);
				VF *backward = model.ComputeBackwardMatrix(seq1, seq2);
				assert(backward);
				// compute posterior probability
				VF *double_posterior = model.ComputePosteriorMatrix(seq1, seq2, *forward,*backward);
				assert(double_posterior);
   				delete forward;
				delete backward;
//global pair-HMM
				// compute posterior probability
				VF *global_posterior = ::ComputePostProbs(a, b, seq1->GetString(),seq2->GetString());
				assert(global_posterior);
//local pair-HMM
				// compute forward and backward probabilities
				forward = model.ComputeForwardMatrix(seq1, seq2,false);
				assert(forward);
				backward = model.ComputeBackwardMatrix(seq1, seq2,false);
				assert(backward);
				// compute posterior probability
				posterior = model.ComputePosteriorMatrix(seq1, seq2, *forward,*backward, false);
				assert(posterior);
   				delete forward;
				delete backward;
//combined model
				//merge probalign + local + probcons
				VF::iterator ptr1 = double_posterior->begin();
				VF::iterator ptr2 = global_posterior->begin();
				VF::iterator ptr = posterior->begin();
				for (int i = 0; i <= seq1->GetLength(); i++) {
					for (int j = 0; j <= seq2->GetLength(); j++) {
						float v1 = *ptr1;
						float v2 = *ptr2;
						float v3 = *ptr;
						//*ptr = v1;
						*ptr = sqrt((v1*v1 + v2*v2 + v3*v3)/3);
						//*ptr = sqrt((v2*v2 + v3*v3)/2);
						ptr1++;
						ptr2++;
						ptr++;
					}
				}
                delete double_posterior;
				delete global_posterior;
			}

            assert(posterior);

			// perform the pairwise sequence alignment
			pair<SafeVector<char> *, float> alignment = model.ComputeAlignment(
			seq1->GetLength(), seq2->GetLength(), *posterior);

			//compute expected accuracy
			distances[a][b] = distances[b][a] = 1.0f - alignment.second
					/ min(seq1->GetLength(), seq2->GetLength());

			// compute sparse representations
			sparseMatrices[a][b] = new SparseMatrix(seq1->GetLength(),
			seq2->GetLength(), *posterior);
			sparseMatrices[b][a] = NULL;

			delete alignment.first;
			delete posterior;
#ifndef _OPENMP
		}
#endif
	}

	timeUsed = GetElapsedTime ( startTime );
	double lastUsed = timeUsed;
	//create the guide tree
	this->tree = new MSAClusterTree(this, distances, numSeqs);
	this->tree->create(vpid);

	// perform the consistency transformation the desired number of times
	for (int r = 0; r < numConsistencyReps; r++) {
		SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices =
				DoRelaxation(sequences, sparseMatrices);
		// now replace the old posterior matrices
		for (int i = 0; i < numSeqs; i++) {
			for (int j = 0; j < numSeqs; j++) {
				delete sparseMatrices[i][j];
				sparseMatrices[i][j] = newSparseMatrices[i][j];
			}
		}
	}
	timeUsed = GetElapsedTime ( startTime );

	//compute the final multiple sequence alignment
	MultiSequence *finalAlignment = ComputeFinalAlignment(this->tree, sequences,
			sparseMatrices, model, pid);

	//cerr << "[Main] MSA alignment Finished : " << endl << endl;

#ifdef _OPENMP
	delete [] seqsPairs;
#endif

	// build annotation
	if (enableAnnotation) {
		WriteAnnotation(finalAlignment, sparseMatrices);
	}
	//destroy the guide tree
	if(this->tree!=NULL) delete this->tree;
	this->tree = 0;

	// delete sparse matrices
	for (int a = 0; a < numSeqs - 1; a++) {
		for (int b = a + 1; b < numSeqs; b++) {
			delete sparseMatrices[a][b];
			delete sparseMatrices[b][a];
		}
	}

	return finalAlignment;
}

//non-progressive alignment
MultiSequence* MSA::npdoAlign( MultiSequence *sequences, const int variance_mean ) {
	assert(sequences);

	//get the number of sequences
	const int numSeqs = sequences->GetNumSequences();
	//creat sparseMatrices
    SafeVector<SafeVector<SparseMatrix *> > sparseMatrices(numSeqs,
			SafeVector<SparseMatrix *>(numSeqs, NULL));
    //create distance matrix
 	VVF distances(numSeqs, VF(numSeqs, 0));

 	startTime = GetTime();
#ifdef _OPENMP
	//calculate sequence pairs for openmp model
	int pairIdx = 0;
	numPairs = (numSeqs - 1) * numSeqs / 2;
	seqsPairs = new SeqsPair[numPairs];
	for(int a = 0; a < numSeqs; a++) {
		for(int b = a + 1; b < numSeqs; b++) {
			seqsPairs[pairIdx].seq1 = a;
			seqsPairs[pairIdx].seq2 = b;
			pairIdx++;
		}
	}
#endif

 	ProbabilisticModel model(initDistrib, gapOpen, gapExtend, emitPairs,emitSingle);

	//compute the posterior pairwise residue alignment probabilities
	ArrangePosteriorProbs(sequences, model, sparseMatrices, distances, variance_mean);

	timeUsed = GetElapsedTime ( startTime );
	double lastUsed = timeUsed;

	// perform the consistency transformation the desired number of times
	for (int r = 0; r < numConsistencyReps; r++) {
		SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices =
				DoRelaxation(sequences, sparseMatrices);
		// now replace the old posterior matrices
		for (int i = 0; i < numSeqs; i++) {
			for (int j = 0; j < numSeqs; j++) {
				delete sparseMatrices[i][j];
				sparseMatrices[i][j] = newSparseMatrices[i][j];
			}
		}
	}

	timeUsed = GetElapsedTime ( startTime );
 	lastUsed = timeUsed;

    MultiSequence *finalAlignment = ComputeGraph(sequences, model, sparseMatrices, distances);
    timeUsed = GetElapsedTime ( startTime );
	lastUsed = timeUsed;

	DoRefinement(finalAlignment, sparseMatrices, model, distances);
	timeUsed = GetElapsedTime ( startTime );

#ifdef _OPENMP
	delete [] seqsPairs;
#endif

	// build annotation
	if (enableAnnotation) {
		WriteAnnotation(finalAlignment, sparseMatrices);
	}
	//destroy the alignment graph
	if(this->graph!=NULL) delete this->graph;
	this->graph=0;

	// delete sparse matrices
	for (int a = 0; a < numSeqs - 1; a++) {
		for (int b = a + 1; b < numSeqs; b++) {
			delete sparseMatrices[a][b];
			delete sparseMatrices[b][a];
		}
	}

	return finalAlignment;
}


/////////////////////////////////////////////////////////////////
// DoRelaxation()
//
// Performs one round of the weighted probabilistic consistency transformation.
//
/////////////////////////////////////////////////////////////////
//no weight
SafeVector<SafeVector<SparseMatrix *> > MSA::DoRelaxation(	MultiSequence *sequences,
		SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices) {
	const int numSeqs = sequences->GetNumSequences();

	SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices(numSeqs,
			SafeVector<SparseMatrix *>(numSeqs, NULL));

	// for every pair of sequences
#ifdef _OPENMP
	int pairIdx;
#pragma omp parallel for private(pairIdx) default(shared) schedule(dynamic)
	for(pairIdx = 0; pairIdx < numPairs; pairIdx++) {
		int i = seqsPairs[pairIdx].seq1;
		int j = seqsPairs[pairIdx].seq2;

#else
	for (int i = 0; i < numSeqs; i++) {
		for (int j = i + 1; j < numSeqs; j++) {
#endif
			Sequence *seq1 = sequences->GetSequence(i);
			Sequence *seq2 = sequences->GetSequence(j);

			if (enableVerbose) {
#ifdef _OPENMP
#pragma omp critical
#endif
			cerr << "Relaxing (" << i + 1 << ") " << seq1->GetHeader()
						<< " vs. " << "(" << j + 1 << ") " << seq2->GetHeader()
						<< ": ";
			}
			// get the original posterior matrix
			VF *posteriorPtr = sparseMatrices[i][j]->GetPosterior();
			assert(posteriorPtr);
			VF &posterior = *posteriorPtr;

			const int seq1Length = seq1->GetLength();
			const int seq2Length = seq2->GetLength();

			// contribution from the summation where z = x and z = y
			for (int k = 0; k < (seq1Length + 1) * (seq2Length + 1); k++) {
				posterior[k] += posterior[k];
			}

			if (enableVerbose)
				cerr << sparseMatrices[i][j]->GetNumCells() << " --> ";

			// contribution from all other sequences
			for (int k = 0; k < numSeqs; k++) {
				if (k != i && k != j) {
					if (k < i)
						Relax1(sparseMatrices[k][i], sparseMatrices[k][j],posterior);
					else if (k > i && k < j)
						Relax(sparseMatrices[i][k], sparseMatrices[k][j],posterior);
					else {
						SparseMatrix *temp = sparseMatrices[j][k]->ComputeTranspose();
						Relax(sparseMatrices[i][k], temp, posterior);
						delete temp;
					}
				}
			}

			for (int k = 0; k < (seq1Length + 1) * (seq2Length + 1); k++) {
				posterior[k] /= numSeqs;
			}
			// mask out positions not originally in the posterior matrix
			SparseMatrix *matXY = sparseMatrices[i][j];
			for (int y = 0; y <= seq2Length; y++)
				posterior[y] = 0;
			for (int x = 1; x <= seq1Length; x++) {
				SafeVector<PIF>::iterator XYptr = matXY->GetRowPtr(x);
				SafeVector<PIF>::iterator XYend = XYptr + matXY->GetRowSize(x);
				VF::iterator base = posterior.begin() + x * (seq2Length + 1);
				int curr = 0;
				while (XYptr != XYend) {
					// zero out all cells until the first filled column
					while (curr < XYptr->first) {
						base[curr] = 0;
						curr++;
					}
					// now, skip over this column
					curr++;
					++XYptr;
				}

				// zero out cells after last column
				while (curr <= seq2Length) {
					base[curr] = 0;
					curr++;
				}
			}

			// save the new posterior matrix
			newSparseMatrices[i][j] = new SparseMatrix(seq1->GetLength(),
					seq2->GetLength(), posterior);
			newSparseMatrices[j][i] = NULL;

			if (enableVerbose)
				cerr << newSparseMatrices[i][j]->GetNumCells() << " -- ";

			delete posteriorPtr;

			if (enableVerbose)
				cerr << "done." << endl;
#ifndef _OPENMP
		}
#endif
	}

	return newSparseMatrices;
}

/////////////////////////////////////////////////////////////////
// Relax()
//
// Computes the consistency transformation for a single sequence
// z, and adds the transformed matrix to "posterior".
/////////////////////////////////////////////////////////////////
//without weight
void MSA::Relax(SparseMatrix *matXZ, SparseMatrix *matZY,
		VF &posterior) {

	assert(matXZ);
	assert(matZY);

	int lengthX = matXZ->GetSeq1Length();
	int lengthY = matZY->GetSeq2Length();
	assert(matXZ->GetSeq2Length() == matZY->GetSeq1Length());

	// for every x[i]
	for (int i = 1; i <= lengthX; i++) {
		SafeVector<PIF>::iterator XZptr = matXZ->GetRowPtr(i);
		SafeVector<PIF>::iterator XZend = XZptr + matXZ->GetRowSize(i);

		VF::iterator base = posterior.begin() + i * (lengthY + 1);

		// iterate through all x[i]-z[k]
		while (XZptr != XZend) {
			SafeVector<PIF>::iterator ZYptr = matZY->GetRowPtr(XZptr->first);
			SafeVector<PIF>::iterator ZYend = ZYptr
					+ matZY->GetRowSize(XZptr->first);
			const float XZval = XZptr->second;

			// iterate through all z[k]-y[j]
			while (ZYptr != ZYend) {
                                base[ZYptr->first] += XZval * ZYptr->second;
				ZYptr++;
			}
			XZptr++;
		}
	}
}

/////////////////////////////////////////////////////////////////
// Relax1()
//
// Computes the consistency transformation for a single sequence
// z, and adds the transformed matrix to "posterior".
/////////////////////////////////////////////////////////////////
//without weight
void MSA::Relax1(SparseMatrix *matZX, SparseMatrix *matZY,VF &posterior) {

	assert(matZX);
	assert(matZY);

	int lengthZ = matZX->GetSeq1Length();
	int lengthY = matZY->GetSeq2Length();

	// for every z[k]
	for (int k = 1; k <= lengthZ; k++) {
		SafeVector<PIF>::iterator ZXptr = matZX->GetRowPtr(k);
		SafeVector<PIF>::iterator ZXend = ZXptr + matZX->GetRowSize(k);

		// iterate through all z[k]-x[i]
		while (ZXptr != ZXend) {
			SafeVector<PIF>::iterator ZYptr = matZY->GetRowPtr(k);
			SafeVector<PIF>::iterator ZYend = ZYptr + matZY->GetRowSize(k);
			const float ZXval = ZXptr->second;
			VF::iterator base = posterior.begin()
					+ ZXptr->first * (lengthY + 1);

			// iterate through all z[k]-y[j]
			while (ZYptr != ZYend) {
				base[ZYptr->first] += ZXval * ZYptr->second;
				ZYptr++;
			}
			ZXptr++;
		}
	}
}

/////////////////////////////////////////////////////////////////
// ProcessTree()
//
// Process the tree recursively.  Returns the aligned sequences
// corresponding to a node or leaf of the tree.
/////////////////////////////////////////////////////////////////

MultiSequence* MSA::ProcessTree(TreeNode *tree, MultiSequence *sequences,
		const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
		const ProbabilisticModel &model) {

	MultiSequence *result;
	// check if this is a node of the alignment tree
	//if (tree->GetSequenceLabel() == -1){
	if (tree->leaf == NODE) {
		MultiSequence *alignLeft = ProcessTree(tree->left, sequences,
				sparseMatrices, model);
		MultiSequence *alignRight = ProcessTree(tree->right, sequences,
				sparseMatrices, model);

		assert(alignLeft);
		assert(alignRight);

		float alignscore;
		result = AlignAlignments(alignLeft, alignRight, sparseMatrices, model, alignscore, true);
		assert(result);

		delete alignLeft;
		delete alignRight;
	}

	// otherwise, this is a leaf of the alignment tree
	else {
		result = new MultiSequence();
		assert(result);
		//result->AddSequence (sequences->GetSequence(tree->GetSequenceLabel())->Clone());
		result->AddSequence(sequences->GetSequence(tree->idx)->Clone());
	}

	return result;
}

/////////////////////////////////////////////////////////////////
// AlignAlignments()
//
// Returns the alignment of two MultiSequence objects.
/////////////////////////////////////////////////////////////////

MultiSequence* MSA::AlignAlignments(MultiSequence *align1,
		MultiSequence *align2,
		const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
		const ProbabilisticModel &model, float &alignscore, bool nflag) {

	// print some info about the alignment
	if (enableVerbose) {
		for (int i = 0; i < align1->GetNumSequences(); i++)
			cerr << ((i == 0) ? "[" : ",")
					<< align1->GetSequence(i)->GetLabel();
		cerr << "] vs. ";
		for (int i = 0; i < align2->GetNumSequences(); i++)
			cerr << ((i == 0) ? "[" : ",")
					<< align2->GetSequence(i)->GetLabel();
		cerr << "]: ";
	}

	VF *posterior;
	if(!nflag)
		posterior = model.BuildPosterior (align1, align2, sparseMatrices, cutoff);
	else
		posterior = model.BuildPosterior(getSeqsWeights(), align1, align2,
			sparseMatrices, cutoff);

	// compute an "accuracy" measure for the MSA before refinement

	pair<SafeVector<char> *, float> alignment;
	//perform alignment
	alignment = model.ComputeAlignment(align1->GetSequence(0)->GetLength(),
			align2->GetSequence(0)->GetLength(), *posterior);

	delete posterior;

	if (enableVerbose) {

		// compute total length of sequences
		int totLength = 0;
		for (int i = 0; i < align1->GetNumSequences(); i++)
			for (int j = 0; j < align2->GetNumSequences(); j++)
				totLength += min(align1->GetSequence(i)->GetLength(),
						align2->GetSequence(j)->GetLength());

		// give an "accuracy" measure for the alignment
		cerr << alignment.second / totLength << endl;
	}

	// now build final alignment
	MultiSequence *result = new MultiSequence();
	for (int i = 0; i < align1->GetNumSequences(); i++)
		result->AddSequence(
				align1->GetSequence(i)->AddGaps(alignment.first, 'X'));
	for (int i = 0; i < align2->GetNumSequences(); i++)
		result->AddSequence(
				align2->GetSequence(i)->AddGaps(alignment.first, 'Y'));
	if (!enableAlignOrder)
		result->SortByLabel();

	// free temporary alignment
	delete alignment.first;
	alignscore = alignment.second;
	return result;
}


/////////////////////////////////////////////////////////////////
// ComputeFinalAlignment()
//
// Compute the final alignment by calling ProcessTree(), then
// performing iterative refinement as needed.
/////////////////////////////////////////////////////////////////

MultiSequence* MSA::ComputeFinalAlignment(MSAGuideTree*tree,
		MultiSequence *sequences,
		const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
		const ProbabilisticModel &model, const int pid) {

	startTime = GetTime();
	MultiSequence *alignment = ProcessTree(tree->getRoot(), sequences,
			sparseMatrices, model);

	SafeVector<int> oldOrdering;
	int numSeqs = alignment->GetNumSequences();
	if (enableAlignOrder) {
		for (int i = 0; i < numSeqs; i++)
			oldOrdering.push_back(alignment->GetSequence(i)->GetSortLabel());
		alignment->SaveOrdering();
		enableAlignOrder = false;
	}

	timeUsed = GetElapsedTime ( startTime );
	double lastUsed = timeUsed;

	if(pid > 3 || numSeqs > 150) numIterativeRefinementReps = 0;
	if(numSeqs <= 50) numIterativeRefinementReps=2*numIterativeRefinementReps;

	int ineffectiveness = 0;
	int i=0;
	int cutoff = 100;

	while( i < numIterativeRefinementReps ) {
		int flag = DoIterativeRefinement(sparseMatrices, model, alignment);

		if( numSeqs > 20 ){
			if(numSeqs < 200){
				if(flag > 0){
			    	if(numIterativeRefinementReps < 4*numSeqs)
						numIterativeRefinementReps ++;
                	if(flag == 1) ineffectiveness ++;
				}

 				//else ineffectiveness = 0;
 				if( ineffectiveness > 2*numSeqs && i>cutoff ) break;
 			}
 			else{
 				if(numSeqs > 200) numIterativeRefinementReps=10;
 			}
		}

		i++;
	}

	timeUsed = GetElapsedTime ( startTime );

	return alignment;
}

//DoIterativeRefinement
int MSA::DoIterativeRefinement(
		const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
		const ProbabilisticModel &model, MultiSequence* &alignment) {
	set<int> groupOne, groupTwo;
	int numSeqs = alignment->GetNumSequences();
	int i;
	// create two separate groups
	for (i = 0; i < numSeqs; i++) {
        int index = rand();
		if (index % 2) {
			groupOne.insert(i);
		} else {
			groupTwo.insert(i);
		}
	}
	if (groupOne.empty() || groupTwo.empty()) return 2;

	// project into the two groups
	MultiSequence *groupOneSeqs = alignment->Project(groupOne);
	assert(groupOneSeqs);
	MultiSequence *groupTwoSeqs = alignment->Project(groupTwo);
	assert(groupTwoSeqs);

//no weight profile-profile for refinement
#if 1
	VF *posterior = model.BuildPosterior (groupOneSeqs, groupTwoSeqs, sparseMatrices, cutoff);
#else
	VF *posterior = model.BuildPosterior(getSeqsWeights(), groupOneSeqs, groupTwoSeqs,
			sparseMatrices, cutoff);
#endif
	// compute an "accuracy" for the currrent MSA before refinement
    SafeVector<SafeVector<char>::iterator> oldOnePtrs(groupOne.size());
	SafeVector<SafeVector<char>::iterator> oldTwoPtrs(groupTwo.size());
    i=0;
	for (set<int>::const_iterator iter = groupOne.begin();
		iter != groupOne.end(); ++iter) {
		oldOnePtrs[i++] = alignment->GetSequence(*iter)->GetDataPtr();
	}

    i=0;
	for (set<int>::const_iterator iter = groupTwo.begin();
		iter != groupTwo.end(); ++iter) {
		oldTwoPtrs[i++] = alignment->GetSequence(*iter)->GetDataPtr();
	}

    VF &posteriorArr = *posterior;
    int oldLength = alignment->GetSequence(0)->GetLength();
	int groupOneindex=0; int groupTwoindex=0;
	float accuracy_before = 0;
    int j;
	for (i = 1; i <= oldLength; i++) {
		// check to see if there is a gap in every sequence of the set
		bool foundOne = false;
		for (j = 0; !foundOne && j < (int) groupOne.size(); j++)
			foundOne = (oldOnePtrs[j][i] != '-');
		// if not, then this column counts towards the sequence length
		if (foundOne) groupOneindex ++;
		bool foundTwo = false;
		for (j = 0; !foundTwo && j < (int) groupTwo.size(); j++)
			foundTwo = (oldTwoPtrs[j][i] != '-');
		if (foundTwo) groupTwoindex ++;
        if(foundOne && foundTwo) accuracy_before +=
				posteriorArr[groupOneindex * (groupTwoSeqs->GetSequence(0)->GetLength() + 1) + groupTwoindex];
	}

	pair<SafeVector<char> *, float> refinealignment;
	//perform alignment
	refinealignment = model.ComputeAlignment(groupOneSeqs->GetSequence(0)->GetLength(),
			groupTwoSeqs->GetSequence(0)->GetLength(), *posterior);
    delete posterior;
	// now build final alignment
	MultiSequence *result = new MultiSequence();
	for (int i = 0; i < groupOneSeqs->GetNumSequences(); i++)
		result->AddSequence(
			groupOneSeqs->GetSequence(i)->AddGaps(refinealignment.first, 'X'));
	for (int i = 0; i < groupTwoSeqs->GetNumSequences(); i++)
		result->AddSequence(
			groupTwoSeqs->GetSequence(i)->AddGaps(refinealignment.first, 'Y'));
	// free temporary alignment
	delete refinealignment.first;
	delete alignment;
    alignment = result;
	delete groupOneSeqs;
	delete groupTwoSeqs;
    if(accuracy_before == refinealignment.second) return 1;
    else return 0;
}


//non-progressive functions
/////////////////////////////////////////////////////////////////
// ArrangePosteriorProbs()
//
// compute the posterior probabilities
// Divergent combined model
// Medium  local model
// Similar  global model
//
/////////////////////////////////////////////////////////////////
void MSA::ArrangePosteriorProbs(MultiSequence *sequences,
		const ProbabilisticModel &model,
		SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
		VVF &distances, const int variance_mean){

	assert(sequences);
	//get the number of sequences
	int pid = variance_mean % 10;

#ifdef _OPENMP
	//calculate sequence pairs for openmp model
	int pairIdx = 0;
#endif

	// do all pairwise alignments for posterior probability matrices
#ifdef _OPENMP
#pragma omp parallel for private(pairIdx) default(shared) schedule(dynamic)
	for(pairIdx = 0; pairIdx < numPairs; pairIdx++) {
		int a= seqsPairs[pairIdx].seq1;
		int b = seqsPairs[pairIdx].seq2;
		if(enableVerbose) {
#pragma omp critical
			cerr <<"tid "<<omp_get_thread_num()<<" a "<<a<<" b "<<b<<endl;
		}
#else
	const int numSeqs = sequences->GetNumSequences();
	for (int a = 0; a < numSeqs - 1; a++) {
		for (int b = a + 1; b < numSeqs; b++) {
#endif
			Sequence *seq1 = sequences->GetSequence(a);
			Sequence *seq2 = sequences->GetSequence(b);

			//posterior probability matrix
			VF* posterior;
			switch (pid) {

//divergent
		    case 0:
			case 1: {
//global pair-HMM
				// compute posterior probability
				VF *global_posterior = ::ComputePostProbs(a, b, seq1->GetString(),seq2->GetString());
				assert(global_posterior);
//local pair-HMM
				// compute forward and backward probabilities
				VF *forward = model.ComputeForwardMatrix(seq1, seq2, false);
				assert(forward);
				VF *backward = model.ComputeBackwardMatrix(seq1, seq2, false);
				assert(backward);
				// compute posterior probability
				posterior = model.ComputePosteriorMatrix(seq1, seq2, *forward,*backward, false);
				assert(posterior);

//double affine pair-HMM
				forward = model.ComputeForwardMatrix(seq1, seq2);
				assert(forward);
				backward = model.ComputeBackwardMatrix(seq1, seq2);
				assert(backward);
				// compute posterior probability
				VF *double_posterior = model.ComputePosteriorMatrix(seq1, seq2, *forward,*backward);
				assert(double_posterior);
//combined model
				//merge probalign + local + probcons
				VF::iterator ptr0 = global_posterior->begin();
				VF::iterator ptr1 = posterior->begin();
				VF::iterator ptr2 = double_posterior->begin();
				for (int i = 0; i <= seq1->GetLength(); i++) {
					for (int j = 0; j <= seq2->GetLength(); j++) {
						float v0 = *ptr0;
						float v1 = *ptr1;
						float v2 = *ptr2;
						//*ptr = v1;
						*ptr1 = sqrt((v0*v0 + v1*v1 + v2*v2)/3);
						//*ptr1 = sqrt((v0*v0 + v1*v1)/2);
						ptr0++;
						ptr1++;
						ptr2++;
					}
				}
				delete forward;
				delete backward;
				delete global_posterior;
				delete double_posterior;
            	assert(posterior);

			}
			 break;
//medium similarity use local pair-HMM
			case 2: {

				VF *forward = model.ComputeForwardMatrix(seq1, seq2 , false);
				assert (forward);
				VF *backward = model.ComputeBackwardMatrix(seq1, seq2 , false);
				assert (backward);
				posterior = model.ComputePosteriorMatrix(seq1, seq2, *forward, *backward, false);
				assert (posterior);
				delete forward;
				delete backward;
				break;
			}

//high similarity use global pair-HMM
			default:
				posterior = ::ComputePostProbs(a, b, seq1->GetString(), seq2->GetString());
				assert(posterior);
				break;
			}

			pair<SafeVector<char> *, float> alignment =
				model.ComputeAlignment(seq1->GetLength(),seq2->GetLength(), *posterior);
			int alignlength = 0;
			for (SafeVector<char>::iterator iter = alignment.first->begin(); iter
					!= alignment.first->end(); ++iter)
				if (*iter == 'B')
					alignlength++;

			float distance = alignment.second / alignlength;//
			distances[a][b] = distances[b][a] = distance;

			// compute sparse representations
			sparseMatrices[a][b] = new SparseMatrix(seq1->GetLength(),
			seq2->GetLength(), *posterior);
			sparseMatrices[b][a] = NULL;

			delete posterior;
#ifndef _OPENMP
		}
#endif
	}
}


/////////////////////////////////////////////////////////////////
// CompuTreeGraph()
//
//   (medium)similar progressive alignment construct guide tree
//	 divergent non-progessive alignment construct alignment graph
//
/////////////////////////////////////////////////////////////////

MultiSequence* MSA:: ComputeGraph(MultiSequence *sequences,
		const ProbabilisticModel &model,
		SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices, VVF &distances) {

	assert(sequences);
	//get the number of sequences


#ifdef _OPENMP
	//calculate sequence pairs for openmp model
	int pairIdx = 0;
#endif

	//Put the pairwise posterior probabilities in "alignp"
	//and the pair residue positions in "aligns"
	VVVI aligns_all(numPairs);
	VVF alignp_all(numPairs);

#ifdef _OPENMP
#pragma omp parallel for private(pairIdx) default(shared) schedule(dynamic)
	for(pairIdx = 0; pairIdx < numPairs; pairIdx++) {
		int a= seqsPairs[pairIdx].seq1;
		int b = seqsPairs[pairIdx].seq2;
		if(enableVerbose) {
#pragma omp critical
			cerr <<"tid "<<omp_get_thread_num()<<" a "<<a<<" b "<<b<<endl;
		}
#else
	int pairIdx = -1;
	const int numSeqs = sequences->GetNumSequences();
	for (int a = 0; a < numSeqs - 1; a++) {
		for (int b = a + 1; b < numSeqs; b++) {
			pairIdx++;
#endif

			SparseMatrix *currSpMat = sparseMatrices[a][b];
			int numRows = currSpMat->GetSeq1Length();

			for (int k = 1; k <= numRows; k++) {
				for (int h = 0; h < currSpMat->GetRowSize(k); h++) {
					VI newentry;
					newentry.push_back(a);
					newentry.push_back(k - 1);
					newentry.push_back(b);
					newentry.push_back((currSpMat->GetRowPtr(k)[h].first) - 1);
					aligns_all[pairIdx].push_back(newentry);
					alignp_all[pairIdx].push_back(currSpMat->GetRowPtr(k)[h].second);
				}
			}
#ifndef _OPENMP
		}
#endif
	}

	VVI aligns;
	VF alignp;
	for(int i= 0; i<numPairs; i++){
		int count = aligns_all[i].size();
		for(int j=0; j<count; j++){
			aligns.push_back(aligns_all[i][j]);
			alignp.push_back(alignp_all[i][j]);
		}
	}

	this->graph = new AlignGraph(aligns, alignp, sequences);
	this->graph->Graph2Align();
	MultiSequence *alignment = this->graph->GetAlignment();
	return alignment;
}

//////////////////////////////////////////////////////////////////
// DoRefinement()
//
// Performs the refinement step
//////////////////////////////////////////////////////////////////

void MSA::DoRefinement(MultiSequence* &alignment, SafeVector<SafeVector<
SparseMatrix *> > &sparseMatrices, const ProbabilisticModel &model, VVF distances) {

/*

	//if(pid == 3) numIterativeRefinementReps=0;
	int numSeqs = alignment->GetNumSequences();
	int ineffectiveness = 0;
	int i=0;
	while( i < numIterativeRefinementReps ) {
		int flag = DoIterativeRefinement(sparseMatrices, model, alignment);
		if( numSeqs > 25 ){
			if(numSeqs < 200){
				if(flag > 0){
			    	if(numIterativeRefinementReps < 4*numSeqs)
						numIterativeRefinementReps ++;
                	if(flag == 1) ineffectiveness ++;
				}

 				//else ineffectiveness = 0;
 				if( ineffectiveness > 2*numSeqs && i>100 ) break;
 			}
 			else{
 				if(numSeqs > 200) numIterativeRefinementReps=10;
 			}
		}
		i++;
	}

*/

	int NumSeqs = alignment->GetNumSequences();
	if(NumSeqs > 150) numIterativeRefinementReps = 0;

	vector<set<int> > SimSeqs;

	FindSimilar(distances, SimSeqs);
	int cnt = 0;
	float oalignscore = 0;
	float nalignscore = 0;
	int ineffectiveness = 0;
	//if(NumSeqs <= 100 ) numIterativeRefinementReps = 2*numIterativeRefinementReps;
	// Performe refinement for at least numRefinementsReAligns number of realignments
	while (cnt < numIterativeRefinementReps) {
		srand(time(0));

		VI list(NumSeqs);

		for (int i = 0; i < NumSeqs; i++)
			list[i] = i;
		VI rnd_list;

		// obtain a random ordering of sequences
		while (list.size()> 0) {
			int index = rand() % (list.size());
			rnd_list.push_back(list[index]);
			list.erase(list.begin() + index);
		}

		//For each sequence, update the set of similar sequences (S_x) and then
		//align that to set of dissimilar sequences (N_x)
		for (int i = 0; i < alignment->GetNumSequences(); i++) {
			int si = rnd_list[i];
			set<int> groupOne, groupTwo;
			groupOne = SimSeqs[si];

			// project sequences to the two groups (S_x,N_x)
			for (int j = 0; j < alignment->GetNumSequences(); j++)
				if (groupOne.find(j) == groupOne.end())
					groupTwo.insert(j);
			cnt++;
			if ((groupOne.size() != 0) & (groupTwo.size() != 0)) {

				MultiSequence *groupOneSeqs = alignment->Project(groupOne);
				assert (groupOneSeqs);
				MultiSequence *groupTwoSeqs = alignment->Project(groupTwo);
				assert (groupTwoSeqs);

				//Find x
				set<int>::iterator it;
				int cnnt = 0;
				for (it = groupOne.begin(); it != groupOne.end(); ++it) {
					if (*it == si)
						break;
					cnnt++;
				}

				float oalignscore2 = 0;
				float nalignscore2 = 0;
				//Update S_x by aligning x with (S_x - x)
				if (groupOneSeqs->GetNumSequences()> 1) {
					set<int> groupOne2, groupTwo2;
					groupOne2.insert(cnnt);

					for (int k = 0; k < groupOneSeqs->GetNumSequences(); k++)
						if (k != cnnt)
							groupTwo2.insert(k);

					MultiSequence *groupOneSeqs2 =
							groupOneSeqs->Project(groupOne2);
					assert (groupOneSeqs2);
					MultiSequence *groupTwoSeqs2 =
							groupOneSeqs->Project(groupTwo2);
					assert (groupTwoSeqs2);
					delete groupOneSeqs;
					// realign
					groupOneSeqs = AlignAlignments(groupOneSeqs2,
							groupTwoSeqs2, sparseMatrices, model, nalignscore2, false);
					if(nalignscore2>oalignscore2) oalignscore2 = nalignscore2;
					else ineffectiveness++;

					cnt++;
				}
				delete alignment;

				// realign the updated similar set (S'_x) and N_x
				alignment = AlignAlignments(groupOneSeqs, groupTwoSeqs,
						sparseMatrices, model, nalignscore, false);
				if(nalignscore<oalignscore && numIterativeRefinementReps<8*NumSeqs && ineffectiveness<4*NumSeqs){
					oalignscore = nalignscore;
					numIterativeRefinementReps += NumSeqs;
				}
			}
		}
	}

}

//////////////////////////////////////////////////////////////////
// FindSimilar()
//
// Find similar sequences for each sequence using kmeans algorithm
//////////////////////////////////////////////////////////////////

void MSA::FindSimilar(VVF distances, vector<set<int> > &SimSeqs) {

	int NumSeqs = distances.size();
	for (int i = 0; i < NumSeqs; i++)
		distances[i][i] = 1;

	// For each sequence x find the set of similar sequences S_x
	for (int i = 0; i < NumSeqs; i++) {
		set<int> c1, c2;

		float min_d = 1;
		float max_d = 0;
		int ii_min = 0;
		int ii_max = 0;
		for (int j = 0; j < NumSeqs; j++) {
			if ((distances[i][j] <= min_d)) {
				ii_min = j;
				min_d = distances[i][j];
			}
			if ((distances[i][j] >= max_d)) {
				ii_max = j;
				max_d = distances[i][j];
			}
		}

		c1.insert(ii_max); // The first cluster: similar sequences (S_x)
		c2.insert(ii_min); // The second cluster: dissimilar sequences (N_x)

		//initiate the clusters
		for (int j = 0; j < NumSeqs; j++) {
			if ((j != ii_min) & (j != ii_max)) {
				if (abs(distances[j][i] - max_d) < abs(distances[j][i] - min_d))
					c1.insert(j);
				else
					c2.insert(j);
			}
		}

		if (c1.find(i) == c1.end()) {
			c2.erase(i);
			c1.insert(i);
		}
		bool ch_flag = true;
		int cnt = 0;

		// Iterate 100 times to obtain the clusters using kmeans algorithm
		while ((cnt < 100) & (ch_flag)) {
			ch_flag = false;
			VI changes(NumSeqs, 0);

			//compute center of each cluster
			float m1 = 0;
			float m2 = 0;
			set<int>::iterator it;
			for (it = c1.begin(); it != c1.end(); ++it)
				m1 += distances[i][*it];
			for (it = c2.begin(); it != c2.end(); ++it)
				m2 += distances[i][*it];
			m1 /= c1.size();
			m2 /= c2.size();

			//update the clusters
			for (int j = 0; j < NumSeqs; j++) {
				if (j != i) {
					set<int>::iterator it;
					if (c1.find(j) != c1.end()) {
						if (abs(distances[j][i] - m1)
								> abs(distances[j][i] - m2)) {
							changes[j] = 1;
							ch_flag = true;
						}
					} else {
						if (abs(distances[j][i] - m2)
								> abs(distances[j][i] - m1)) {
							changes[j] = -1;
							ch_flag = true;
						}
					}
				}
			}
			if (ch_flag) {
				for (int j = 0; j < NumSeqs; j++) {
					if (changes[j] == 1) {
						c1.erase(j);
						c2.insert(j);
					} else if (changes[j] == -1) {
						c2.erase(j);
						c1.insert(j);
					}
				}
			}
			cnt++;
		}
		SimSeqs.push_back(c1);
	}

}



/////////////////////////////////////////////////////////////////
// GetInteger()
//
// Attempts to parse an integer from the character string given.
// Returns true only if no parsing error occurs.
/////////////////////////////////////////////////////////////////

bool MSA::GetInteger(char *data, int *val) {
	char *endPtr;
	long int retVal;

	assert(val);

	errno = 0;
	retVal = strtol(data, &endPtr, 0);
	if (retVal == 0 && (errno != 0 || data == endPtr))
		return false;
	if (errno != 0 && (retVal == LONG_MAX || retVal == LONG_MIN))
		return false;
	if (retVal < (long) INT_MIN || retVal > (long) INT_MAX)
		return false;
	*val = (int) retVal;
	return true;
}

/////////////////////////////////////////////////////////////////
// GetFloat()
//
// Attempts to parse a float from the character string given.
// Returns true only if no parsing error occurs.
/////////////////////////////////////////////////////////////////

bool MSA::GetFloat(char *data, float *val) {
	char *endPtr;
	double retVal;

	assert(val);

	errno = 0;
	retVal = strtod(data, &endPtr);
	if (retVal == 0 && (errno != 0 || data == endPtr))
		return false;
	if (errno != 0 && (retVal >= 1000000.0 || retVal <= -1000000.0))
		return false;
	*val = (float) retVal;
	return true;
}


/////////////////////////////////////////////////////////////////
// WriteAnnotation()
//
// Computes annotation for multiple alignment and write values
// to a file.
/////////////////////////////////////////////////////////////////

void MSA::WriteAnnotation(MultiSequence *alignment,
		const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices) {
	ofstream outfile(annotationFilename.c_str());

	if (outfile.fail()) {
		cerr << "ERROR: Unable to write annotation file." << endl;
		exit(1);
	}

	const int alignLength = alignment->GetSequence(0)->GetLength();
	const int numSeqs = alignment->GetNumSequences();

	SafeVector<int> position(numSeqs, 0);
	SafeVector<SafeVector<char>::iterator> seqs(numSeqs);
	for (int i = 0; i < numSeqs; i++)
		seqs[i] = alignment->GetSequence(i)->GetDataPtr();
	SafeVector<pair<int, int> > active;
	active.reserve(numSeqs);

	SafeVector<int> lab;
	for (int i = 0; i < numSeqs; i++)
		lab.push_back(alignment->GetSequence(i)->GetSortLabel());

	// for every column
	for (int i = 1; i <= alignLength; i++) {

		// find all aligned residues in this particular column
		active.clear();
		for (int j = 0; j < numSeqs; j++) {
			if (seqs[j][i] != '-') {
				active.push_back(make_pair(lab[j], ++position[j]));
			}
		}

		sort(active.begin(), active.end());
		outfile << setw(4) << ComputeScore(active, sparseMatrices) << endl;
	}

	outfile.close();
}

/////////////////////////////////////////////////////////////////
// ComputeScore()
//
// Computes the annotation score for a particular column.
/////////////////////////////////////////////////////////////////

int MSA::ComputeScore(const SafeVector<pair<int, int> > &active,
		const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices) {

	if (active.size() <= 1)
		return 0;

	// ALTERNATIVE #1: Compute the average alignment score.

	float val = 0;
	for (int i = 0; i < (int) active.size(); i++) {
		for (int j = i + 1; j < (int) active.size(); j++) {
			val += sparseMatrices[active[i].first][active[j].first]->GetValue(
					active[i].second, active[j].second);
		}
	}

	return (int) (200 * val / ((int) active.size() * ((int) active.size() - 1)));

}
