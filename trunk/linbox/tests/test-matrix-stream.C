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

/* 
I would like to see a matrix writer that 
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
bool testMatrixStream(std::ostream& out) {
        std::vector<const char *> matrixNames;
	matrixNames. push_back("data/sms.matrix");
	matrixNames.push_back("data/matrix-market-array.matrix");
	matrixNames.push_back("data/maple-sparse1.matrix");
	matrixNames.push_back("data/maple-dense1.matrix");
	matrixNames.push_back("data/generic-dense.matrix");
	matrixNames.push_back("data/sparse-row.matrix");
	matrixNames.push_back("data/matrix-market-coordinate.matrix");

	bool fail = false;
	bool failThis;
	commentator.start("Testing matrix-stream...", "matrix stream", 1);

	for (size_t i = 0; i < 19; ++i) {matrix[2][2] *= 10; matrix[2][2] += 8; }

	for( size_t i = 0; i < matrixNames.size(); ++i ) {
		failThis = false;
		std::ifstream fin(matrixNames[i]);
		if( !fin ) {
			out << "Could not open " << matrixNames[i] << std::endl;
			fail = true;
			continue;
		}
		out << "\tTesting " << matrixNames[i] << std::endl;
		MatrixStream<TestField > ms(f,fin);
		int nzCount = nonZeros;
		size_t m, n;
		integer v;
		if(!ms.getDimensions(m,n)) {
			out << "Error getting dimensions in "
			    << matrixNames[i] << std::endl;
			fail = failThis = true;
		}
		if( !failThis && m != rowDim ) {
			out << "Wrong rowDim in " << matrixNames[i]
			     << ", format " << ms.getFormat() << std::endl
			     << "Got " << m << ", should be " << rowDim << std::endl;
			fail = failThis = true;
		}
		if( !failThis && n != colDim ) {
			out << "Wrong colDim in " << matrixNames[i]
			     << ", format " << ms.getFormat() << std::endl
			     << "Got " << n << ", should be " << colDim << std::endl;
			fail = failThis = true;
		}
		while( !failThis && ms.nextTriple(m,n,v) ) {
			if(!f.isZero(v)) --nzCount;
			if( m >= rowDim ) {
				out << "Row index out of bounds in "
				     << matrixNames[i]
				     << ", format " << ms.getFormat() << std::endl
				     << "Got " << m << ", should be less than "
				     << rowDim << std::endl;
				fail = failThis = true;
				break;
			}
			if( n >= colDim ) {
				out << "Column index out of bounds in "
				     << matrixNames[i]
				     << ", format " << ms.getFormat() << std::endl
				     << "Got " << n << ", should be less than "
				     << colDim << std::endl;
				fail = failThis = true;
				break;
			}
			if( matrix[m][n] != v ) {
				out << "Invalid entry in "
				     << matrixNames[i]
				     << ", format " << ms.getFormat() << std::endl
				     << "Got " << v << " at index (" << m
				     << "," << n << "), should be "
				     << matrix[m][n] << std::endl;
				fail = failThis = true;
				break;
			}
		}
		if( ms.getError() != END_OF_MATRIX ) {
			fail = failThis = true;
		}
		if( !failThis && nzCount > 0 ) {
			out << "Not enough entries in " << matrixNames[i]
			     << ", format " << ms.getFormat() << std::endl;
			fail = failThis = true;
		}
		if( !failThis && nzCount < 0 ) {
			out << "Duplicate entries in " << matrixNames[i]
			     << ", format " << ms.getFormat() << std::endl;
			fail = failThis = true;
		}

		if( failThis ) {
			std::cout << "Test failed for " << matrixNames[i]
			     << ", format " << ms.getFormat() << std::endl
			     << "Error code: " << ms.getError()
			     << ", line number: " << ms.getLineNumber() << std::endl;
			continue;
		}

		fin.seekg(0,std::ios::beg);
		MatrixStream<TestField> ms2(f,fin);
		std::vector<TestField::Element> array;
		if( !ms2.getArray(array) ) {
			ms2.reportError("Problem with getArray",ms2.getLineNumber());
			fail = failThis = true;
		}
		if( !failThis && array.size() != rowDim*colDim ) {
			out << "Array given wrong size in " << matrixNames[i]
			     << ", format " << ms2.getFormat() << std::endl
			     << "Got " << array.size() << ", should be "
			     << rowDim*colDim << std::endl;
			fail = failThis = true;
		}
		for( m = 0; !failThis && m < rowDim; ++m ) {
			for( n = 0; !failThis && n < colDim; ++n ) {
				if( array[m*colDim+n] != matrix[m][n] ) {
					out << "Invalid entry in getArray of "
						  << matrixNames[i]
						  << ", format " << ms2.getFormat()
						  << std::endl
						  << "Got " << array[m*colDim+n] 
						  << " at index (" << m
						  << "," << n << "), should be "
						  << matrix[m][n] << std::endl;
					fail = failThis = true;
					break;
				}
			}
		}

		if( failThis ) {
			out << "Test failed for " << matrixNames[i]
			     << ", format " << ms2.getFormat() << std::endl
			     << "Error code: " << ms2.getError()
			     << ", line number: " << ms2.getLineNumber() << std::endl;
		}
	}

	if( 	testMatrix< DenseMatrix<TestField> >
			( out, matrixNames[0], "Dense BlackBox Matrix" )
	  ) fail = true;
	if( 	testMatrix< SparseMatrix<TestField> >
			( out, matrixNames[0], "Sparse BlackBox Matrix" )
	  ) fail = true;
	if( 	testMatrix< BlasBlackbox<TestField> >
			( out, matrixNames[0], "BLAS BlackBox Matrix" )
	  ) fail = true;
	
	if( fail )	out << "FAIL: matrix-stream" << std::endl; 
	else 		out << "matrix-stream Passed" << std::endl;
	return !fail;
}

template <class BB>
bool testMatrix( std::ostream& out, const char* filename, const char* BBName ) {
	out << "\tTesting " << BBName << std::endl;
	std::ifstream fin( filename );
	if( !fin ) {
		out << "Could not open " << filename << std::endl;
		return true;
	}
	MatrixStream<TestField > ms(f, fin);
	BB m( ms );
	if( m.rowdim() != rowDim ) {
		out << "Wrong rowDim in " << BBName << std::endl
		     << "Got " << m.rowdim() << ", should be " << rowDim << std::endl;
		return true;
	}
	if( m.coldim() != colDim ) {
		out << "Wrong colDim in " << BBName << std::endl
		     << "Got " << m.coldim() << ", should be " << colDim << std::endl;
		return true;
	}
	bool fail = false;
	for( size_t i = 0; i < rowDim; ++i ) {
		for( size_t j = 0; j < colDim; ++j ) {
			if( m.getEntry(i,j) != matrix[i][j] ) {
				out << "Invalid entry in " << BBName << " at index ("
				     << i << "," << j << ")" << std::endl
				     << "Got " << m.getEntry(i,j) << ", should be "
				     << matrix[i][j] << std::endl;
				fail = true;
			}
		}
	}
	return fail;
}

int main(int argc, char* argv[]){
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

	commentator.start("Matrix stream test suite", "matrix stream", 1);
	std::ostream& report = commentator.report();
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	bool pass = testMatrixStream(report);
	commentator.stop("matrix stream test suite");
	return pass ? 0 : -1;
}
