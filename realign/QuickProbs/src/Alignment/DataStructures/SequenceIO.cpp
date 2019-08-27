#include <vector>
#include <stdexcept>
#include <limits>
#include <cstring>
#include <algorithm>

#include "Sequence.h"
#include "MultiSequence.h"
#include "SequenceIO.h"

using namespace std;
using namespace quickprobs;

/// <summary>
/// See declaration for all the details.
/// </summary>
void quickprobs::SequenceIO::load(std::string filename, SequenceFormat format, MultiSequence& set)
{
	std::ifstream file;
	file.open(filename.c_str(), ios_base::binary);

	if (!file.is_open()) {
		throw std::runtime_error("SequenceIO::load(): unable to open input file.");
	}

	switch (format) {
	case FASTA:
		loadFasta(file, set); break;
	case CLUSTALW:
		loadClustal(file, set); break;
	case PHYLIP:
		loadPhylip(file, set); break;
	default:
		throw std::runtime_error("SequenceIO::load(): unrecognised file format.");
	}

	bool ok = checkAndCorrect(set);

	if (!ok) {
		throw std::runtime_error("Illegal characters in sequence set!");
	}
}


void quickprobs::SequenceIO::save(std::string filename, SequenceFormat format, MultiSequence& set)
{
	std::ofstream file;
	file.open(filename.c_str(), ios_base::binary);

	if (!file.is_open()) {
		throw std::runtime_error("SequenceIO::save(): unable to open output file.");
	}

	switch (format) {
	case FASTA:
		saveFasta(file, set); break;
	case CLUSTALW:
		saveClustal(file, set); break;
	case PHYLIP:
		savePhylip(file, set); break;
	default:
		throw std::runtime_error("SequenceIO::save(): unrecognised file format.");
	}
}


/// <summary>
/// See declaration for all the details.
/// </summary>
bool quickprobs::SequenceIO::checkAndCorrect(MultiSequence& set)
{
	bool ok = true;

	for (int i = 0; i < set.count(); ++i) {
		auto seq = set.GetSequence(i);
		char *data = seq->getData() + 1;
		int length = seq->GetLength();

		// control and correct characters (omit @ at the beginning)
		std::transform(data, data + length, data, [&ok](char c)->char {
			if		(c == '.')		{ c = '-'; }
			else if	(isalpha(c))	{ c = toupper(c); }
			else {
				cout << "illegal sequence character:" << c << endl;
				ok = false;
			}

			return c;
		});
	}

	return ok;
}

/// <summary>
/// See declaration for all the details.
/// </summary>
void quickprobs::SequenceIO::loadFasta(std::istream& input, MultiSequence& set)
{
	char buffer[MAX_LINE_LENGTH + 1];
	int added = 0;

	// read all sequences
	while (input.good()) {
		input.getline(buffer, MAX_LINE_LENGTH);

		// omit empty lines
		if (strlen(buffer) == 0) { continue; }

		// check to make sure that it is a correct header line
		if (buffer[0] == '>') {
			std::string header;
			header.assign(buffer + 1, buffer + strlen(buffer)); // remove the leading ">"

			// remove any leading or trailing white space in the header comment
			while (header.length() > 0 && isspace (header[0])) { header = header.substr(1); }
			while (header.length() > 0 && isspace (header[header.length() - 1])) { header = header.substr(0, header.length() - 1); }

			// get the sequence label as being the current # of sequences (sequence labels here are zero-based)
			int index = set.count();
			auto data = new std::vector<char>;
			data->reserve(MAX_LINE_LENGTH);
			data->push_back('@');

			// get consecutive sequence
			while (input.good()) {
				// get the first character;
				if (input.peek() != '>') {
					input.getline(buffer, MAX_LINE_LENGTH);

					// omit empty lines
					if (strlen(buffer) == 0) { continue; }

					// remove \r if exists
					if (buffer[strlen(buffer) - 1] == '\r') { buffer[strlen(buffer) - 1] = 0; }

					data->insert(data->end(), buffer, buffer + strlen(buffer));
				} else {
					break;
				}
			}

			//data->push_back(0);

			Sequence *seq = new Sequence(data, header, data->size() - 1, index, index);
			++added;

			set.AddSequence(seq);
		}
	}

	if (added == 0) {
		throw std::runtime_error("SequenceIO::loadFasta(): no sequences read.");
	}
}

/// <summary>
/// See declaration for all the details.
/// </summary>
void quickprobs::SequenceIO::loadClustal(std::istream& input, MultiSequence& set)
{
	throw std::runtime_error("SequenceIO::loadClustal(): function not implemented.");
}

/// <summary>
/// See declaration for all the details.
/// </summary>
void quickprobs::SequenceIO::loadPhylip(std::istream& input, MultiSequence& set)
{
	throw std::runtime_error("SequenceIO::loadClustal(): function not implemented.");
}

/// <summary>
/// See declaration for all the details.
/// </summary>
void quickprobs::SequenceIO::saveFasta(std::ostream& output, const MultiSequence& set)
{
	// loop through all sequences and write them out
	for (int i = 0; i != set.count(); ++i) {
		auto seq = set.GetSequence(i);
		output << ">" << seq->getHeader() << endl;

		// print out character data
		int ct;
		for (ct = 1; ct <= seq->GetLength(); ct++){
			output << seq->getIterator()[ct];
			if (ct % FASTA_LINE_LENGTH == 0) {
				output << endl;
			}
		}
		if ((ct - 1) % FASTA_LINE_LENGTH != 0) {
			output << endl;
		}
	}
}

/// <summary>
/// See declaration for all the details.
/// </summary>
void quickprobs::SequenceIO::saveClustal(std::ostream& output, const MultiSequence& set)
{
	throw std::runtime_error("SequenceIO::saveClustal(): function not implemented.");
}


/// <summary>
/// See declaration for all the details.
/// </summary>
void SequenceIO::savePhylip(std::ostream& output, const MultiSequence& alignment)
{
	int count = alignment.count();
	int length = alignment.GetSequence(0)->GetLength(); // assume all sequences to be of equal length;

	output << count << endl;
	output << length << endl;

	int globalPosition = 0;

	for (int i = 0; i < count; ++i) {
		auto seq = alignment.GetSequence(i);
		const char* data = seq->getData() + 1;
		std::string header = seq->getHeader();
		header.resize(10, ' ');

		output << header;

		int j;
		for (j = 0; j < std::min(length, 50); ++j) {
			if (j % 10 == 0) {
				output << " ";
			}

			output << data[j];
		}
		output << endl;

		globalPosition = j;
	}
	output << endl;

	while (globalPosition < length) {

		int sequencePosition;
		for (int i = 0; i < count; ++i) {
			auto seq = alignment.GetSequence(i);
			const char* data = seq->getData() + 1;

			output << std::string(10, ' ');

			int k = 0;
			for (sequencePosition = globalPosition; (sequencePosition < length) && (k < 50); ++sequencePosition, ++k) {
				if (sequencePosition % 10 == 0) {
					output << " ";
				}

				output << data[sequencePosition];
			}
			output << endl;
		}

		globalPosition = sequencePosition;
		output << endl;
	}
}



