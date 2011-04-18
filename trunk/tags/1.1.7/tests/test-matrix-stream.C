/* Copyright (C) LinBox
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "test-common.h"
#include <linbox/field/unparametric.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/integer.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/blackbox/blas-blackbox.h>

using namespace LinBox;

/** @file tests/test-matrix-stream.C
@todo I would like to see a matrix writer that 
writes sms format and generic dense format.  
Then we could have a self contained test that 
checks the write/read cycle without depending on preexisting data files.

...but data files illustrating formats that we intend to read but not write would continue to be used.
-bds
*/

const size_t rowDim = 11;
const size_t colDim = 11;
int nonZeros = 33;
integer matrix[rowDim][colDim] = {
				{0, 0, 2, 3, 0, 0, 0, 0, 0, 1, 0},
				{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
				{2, 0, 0, 1, -1, 0, 0, 0, 0, 0, 6},
				{3, 0, 1, 4, 0, 0, 12, 0, 0, -13, 0},
				{0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0},
				{0, 0, 0, 12, 0, 0, 0, 0, 0, 0, 0},
				{0, 0, 0, 0, 0, 1, 0, 500, 400, 300, 200},
				{0, 0, 0, 0, 0, 0, 0, 400, 0, 0, 0},
				{1, 0, 0, -13, 0, 1, 0, 300, 0, 10, 1},
				{0, 0, 6, 0, 0, 0, 0, 200, 0, 1, 0} };
typedef UnparametricField<integer> TestField;
TestField f;

template <class BB>
bool testMatrix( std::ostream& out, const char* filename, const char* BBName ) ;
bool testMatrixStream(const string& matfile) 
{

	bool pass = true;
	commentator.start("Testing matrix-stream...", matfile.c_str());
	std::ostream& out = commentator.report();

	std::ifstream fin(matfile.c_str());
	if( !fin ) {
		out << "Could not open " << matfile << std::endl;
		pass = false;
		}
	out << "\tTesting " << matfile << std::endl;
	MatrixStream<TestField > ms(f,fin);
	int nzCount = nonZeros;
	size_t i, j;
	integer v;
	if(!ms.getDimensions(i,j)) {
		out << "Error getting dimensions in "
		    << matfile << std::endl;
		pass = false;
	}
	if( pass && i != rowDim ) {
		out << "Wrong rowDim in " << matfile
		     << ", format " << ms.getFormat() << std::endl
		     << "Got " << i << ", should be " << rowDim << std::endl;
		pass = false;
	}
	if( pass && j != colDim ) {
		out << "Wrong colDim in " << matfile
		     << ", format " << ms.getFormat() << std::endl
		     << "Got " << j << ", should be " << colDim << std::endl;
		pass = false;
	}
	while( pass && ms.nextTriple(i,j,v) ) {
		if(!f.isZero(v)) --nzCount;
		if( i >= rowDim ) {
			out << "Row index out of bounds in "
			     << matfile
			     << ", format " << ms.getFormat() << std::endl
			     << "Got " << i << ", should be less than "
			     << rowDim << std::endl;
			pass = false;
			break;
		}
		if( j >= colDim ) {
			out << "Column index out of bounds in "
			     << matfile
			     << ", format " << ms.getFormat() << std::endl
			     << "Got " << j << ", should be less than "
			     << colDim << std::endl;
			pass = false;
			break;
		}
		if( matrix[i][j] != v ) {
			out << "Invalid entry in "
			     << matfile
			     << ", format " << ms.getFormat() << std::endl
			     << "Got " << v << " at index (" << i
			     << "," << j << "), should be "
			     << matrix[i][j] << std::endl;
			pass = false;
			break;
		}
	}
	if( ms.getError() != END_OF_MATRIX ) {
		pass = false;
	}
	if( pass && nzCount > 0 ) {
		out << "Not enough entries in " << matfile
		     << ", format " << ms.getFormat() << std::endl;
		pass = false;
	}
	if( pass && nzCount < 0 ) {
		out << "Duplicate entries in " << matfile
		     << ", format " << ms.getFormat() << std::endl;
		pass = false;
	}

	if( !pass ) {
		out << "Test failed for " << matfile
		     << ", format " << ms.getFormat() << std::endl
		     << "Error code: " << ms.getError()
		     << ", line number: " << ms.getLineNumber() << std::endl;
		pass = false;
	}

	fin.seekg(0,std::ios::beg);
	MatrixStream<TestField> ms2(f,fin);
	std::vector<TestField::Element> array;
	if( !ms2.getArray(array) ) {
		ms2.reportError("Problem with getArray",ms2.getLineNumber());
		pass = false;
	}
	if( pass && array.size() != rowDim*colDim ) {
		out << "Array given wrong size in " << matfile
		     << ", format " << ms2.getFormat() << std::endl
		     << "Got " << array.size() << ", should be "
		     << rowDim*colDim << std::endl;
		pass = false;
	}
	for( i = 0; pass && i < rowDim; ++i ) {
		for( j = 0; pass && j < colDim; ++j ) {
			if( array[i*colDim+j] != matrix[i][j] ) {
				out << "Invalid entry in getArray of "
					  << matfile
					  << ", format " << ms2.getFormat()
					  << std::endl
					  << "Got " << array[i*colDim+j] 
					  << " at index (" << i
					  << "," << j << "), should be "
					  << matrix[i][j] << std::endl;
				pass = false;
				break;
			}
		}
	}

	if( !pass ) {
		out << "Test failed for " << matfile
		     << ", format " << ms2.getFormat() << std::endl
		     << "Error code: " << ms2.getError()
		     << ", line number: " << ms2.getLineNumber() << std::endl;
	}

/* later
	if( !testMatrix< DenseMatrix<TestField> >
			( out, matfile[0], "Dense BlackBox Matrix" )
	  ) pass = false;
	if( !testMatrix< SparseMatrix<TestField> >
			( out, matfile[0], "Sparse BlackBox Matrix" )
	  ) pass = false;
	if( !testMatrix< BlasBlackbox<TestField> >
			( out, matfile[0], "BLAS BlackBox Matrix" )
	  ) pass = false;
	
*/
	if( !pass )	out << "FAIL: matrix-stream" << std::endl; 
	else 		out << "matrix-stream Passed" << std::endl;
	commentator.stop(ms.getFormat());
	//commentator.stop(MSG_STATUS(pass));
	return pass;
}

template <class BB>
bool testMatrix( std::ostream& out, const char* filename, const char* BBName ) 
{
	out << "\tTesting " << BBName << std::endl;
	std::ifstream fin( filename );
	if( !fin ) {
		out << "Could not open " << filename << std::endl;
		return false;
	}
	MatrixStream<TestField > ms(f, fin);
	BB m( ms );
	if( m.rowdim() != rowDim ) {
		out << "Wrong rowDim in " << BBName << std::endl
		     << "Got " << m.rowdim() << ", should be " << rowDim << std::endl;
		return false;
	}
	if( m.coldim() != colDim ) {
		out << "Wrong colDim in " << BBName << std::endl
		     << "Got " << m.coldim() << ", should be " << colDim << std::endl;
		return false;
	}
	bool pass = true;
	for( size_t i = 0; i < rowDim; ++i ) {
		for( size_t j = 0; j < colDim; ++j ) {
			if( m.getEntry(i,j) != matrix[i][j] ) {
				out << "Invalid entry in " << BBName << " at index ("
				     << i << "," << j << ")" << std::endl
				     << "Got " << m.getEntry(i,j) << ", should be "
				     << matrix[i][j] << std::endl;
				pass = false;
			}
		}
	}
	return pass;
}

int main(int argc, char* argv[])
{
/*
    static size_t n = 10;
    static integer q = 4093U;
    static int iterations = 2;

*/
    static Argument args[] = {
	/*
    { 'n', "-n N", "Set dimension of test matrices to NxN", TYPE_INT,     &n },
	{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1]", TYPE_INTEGER, &q },
	{ 'i', "-i I", "Perform each test for I iterations",    TYPE_INT,     &iterations },
	{ '\0' }
	*/
	};

	parseArguments (argc, argv, args);

	for (size_t i = 0; i < 19; ++i) {matrix[2][2] *= 10; matrix[2][2] += 8; }

	commentator.start("Matrix stream test suite", "matrix stream");
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	bool pass = true;
	pass = pass && testMatrixStream("data/sms.matrix");
	pass = pass && testMatrixStream("data/matrix-market-array.matrix");
	pass = pass && testMatrixStream("data/maple-sparse1.matrix");
	pass = pass && testMatrixStream("data/maple-dense1.matrix");
	pass = pass && testMatrixStream("data/generic-dense.matrix");
	pass = pass && testMatrixStream("data/sparse-row.matrix");
	pass = pass && testMatrixStream("data/matrix-market-coordinate.matrix");
	commentator.stop(MSG_STATUS(pass));
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
