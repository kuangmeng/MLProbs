#pragma once
#include <limits>
#include <stdexcept>

#include "ISparseMatrix.h"
#include "Common/Printable.h"

class PackedSparseMatrix : public FixedSparseMatrix, public Printable
{
public:
	typedef FixedSparseMatrix BaseType; 

	static int sparseThreshold;
	static void setSparseRowThreshold(int threshold) { sparseThreshold = threshold; }

	const index_type* getRowSizes() const		{ return rowSizes; }
	const index_type* getRowIndices() const		{ return rowIndices; }
	const column_type* getColumns() const		{ throw std::runtime_error("Not implemented!"); return nullptr; }
	const value_type* getValues() const			{ throw std::runtime_error("Not implemented!"); return nullptr; }
	const cell_type* getData() const			{ return data; }
	const cell_type* getRowPtr (int row) const	{ return data + rowIndices[row]; }

	int getRowSize (int row) const				{ return rowSizes[row]; }
	int getRowIndex (int row) const				{ return rowIndices[row]; }
	int getSeq1Length () const					{ return height - 1; }
	int getSeq2Length () const					{ return width - 1; }
	int getHeight() const						{ return height; }
	int getWidth() const						{ return width; }
	int getNumCells () const					{ return numCells; }
	
	int getMaxRow() const						{ return rowSizes[0]; }


	const ::size_t metaBytesNeeded() const { return BaseType::metaBytesNeeded(getSeq1Length()); }
	const ::size_t bytesNeeded() const { return BaseType::bytesNeeded(getSeq1Length(), getNumCells()); }

	float getValue (int row, int col) const;

	PackedSparseMatrix(int seq1Length, int seq2Length, ::size_t numCells);

	PackedSparseMatrix(int seq1Length, int seq2Length, const float*  posterior, float cutoff);

	~PackedSparseMatrix();

	template <class InMatrixType>
	float unpack(const buffer_type* ptr);

	template <class OutMatrixType>
	void pack(buffer_type *ptr, bool freeBuffer, float initialDistance);

	PackedSparseMatrix* computeTranspose() const;

	void fillTransposed(const PackedSparseMatrix& ref);

	std::vector<float> *getPosterior () const;

	void shrinkToFit();

	virtual std::string toString();

protected:

	buffer_type *buffer;

	int width;
	int height;
	int numCells;

	index_type* rowSizes;
	index_type* rowIndices;
	cell_type* data;

	PackedSparseMatrix(int seq1Length, int seq2Length) 
		: height(seq1Length + 1), width(seq2Length + 1), numCells(0), buffer(nullptr) {}

};

#include "PackedSparseMatrix.hpp"