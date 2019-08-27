#include <fstream>
#include <algorithm>
#include <stdexcept>

#include "PhylipTree.h"
#include "DataStructures/MultiSequence.h"
#include "DataStructures/SequenceIO.h"

using namespace std;
using namespace quickprobs;

PhylipTree::PhylipTree(quickprobs::MultiSequence& sequences, PhylogenyType type) 
	: NewickTree(sequences, ""), type(type)
{
	
}

void PhylipTree::build()
{
	std::ofstream seqfile;
	std::ofstream cfgfile;
	std::ostringstream oss;

	seqfile.open("infile", ofstream::out);
	cfgfile.open("config", ofstream::out);
	
	if (!seqfile.is_open() || !cfgfile.is_open()) {
		runtime_error("Unable to create Phylip files.");
	}

	if (type == LIKELIHOOD) {
		cfgfile << "P" << endl << "P" << endl << "Y";
		cfgfile.close();
		oss << "proml.exe < config";
	} else if (type == PARSIMONY) {
		oss << "protpars.exe < config";
	}

	// save alignment as Phylip file
	std::remove("outfile");
	std::remove("outtree");
	SequenceIO::savePhylip(seqfile, sequences);

	// execute Phylip
	std::system(oss.str().c_str());

	// load output tree
	std::ifstream treefile;
	treefile.open("outtree", ifstream::in | ifstream::binary );
	
	std::streampos fsize = 0;
	treefile.seekg( 0, ifstream::end);
	fsize = treefile.tellg();
	treefile.seekg(0, ifstream::beg);
	
	char* buffer = new char[fsize];
	treefile.read(buffer, fsize);
	description.assign(buffer, buffer + fsize);
	description.erase(std::remove(description.begin(), description.end(), 13), description.end());
	description.erase(std::remove(description.begin(), description.end(), 10), description.end());
	delete [] buffer;

	NewickTree::build();
}


