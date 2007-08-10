#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <linbox/field/unparametric.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/integer.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/blackbox/blas-blackbox.h>

using namespace LinBox;

const int nMatrices = 7;
const char* matrixNames[nMatrices] = {
				"data/sms.matrix",
				"data/matrix-market-array.matrix",
				"data/maple-sparse1.matrix",
				"data/maple-dense1.matrix",
				"data/generic-dense.matrix",
				"data/sparse-row.matrix",
				"data/matrix-market-coordinate.matrix"
			       };

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
bool testBlackBox( const char* filename, const char* BBName ) {
	std::cout << "\tTesting " << BBName << std::endl;
	std::ifstream fin( filename );
	if( !fin ) {
		std::cout << "Could not open " << filename << std::endl;
		return true;
	}
	MatrixStream<TestField > ms(f, fin);
	BB m( ms );
	if( m.rowdim() != rowDim ) {
		std::cout << "Wrong rowDim in " << BBName << std::endl
		     << "Got " << m.rowdim() << ", should be " << rowDim << std::endl;
		return true;
	}
	if( m.coldim() != colDim ) {
		std::cout << "Wrong colDim in " << BBName << std::endl
		     << "Got " << m.coldim() << ", should be " << colDim << std::endl;
		return true;
	}
	bool fail = false;
	for( size_t i = 0; i < rowDim; ++i ) {
		for( size_t j = 0; j < colDim; ++j ) {
			if( m.getEntry(i,j) != matrix[i][j] ) {
				std::cout << "Invalid entry in " << BBName << " at index ("
				     << i << "," << j << ")" << std::endl
				     << "Got " << m.getEntry(i,j) << ", should be "
				     << matrix[i][j] << std::endl;
				fail = true;
			}
		}
	}
	return fail;
}

int main() {
	bool fail = false;
	bool failThis;
	std::cout << "Testing matrix-stream..." << std::endl;

	for (int i = 0; i < 19; ++i) {matrix[2][2] *= 10; matrix[2][2] += 8; }

	for( int i = 0; i < nMatrices; ++i ) {
		failThis = false;
		std::ifstream fin(matrixNames[i]);
		if( !fin ) {
			std::cout << "Could not open " << matrixNames[i] << std::endl;
			fail = true;
			continue;
		}
		std::cout << "\tTesting " << matrixNames[i] << std::endl;
		MatrixStream<TestField > ms(f,fin);
		int nzCount = nonZeros;
		size_t m, n;
		integer v;
		if(!ms.getDimensions(m,n)) {
			std::cout << "Error getting dimensions in "
			          << matrixNames[i] << std::endl;
			fail = failThis = true;
		}
		if( !failThis && m != rowDim ) {
			std::cout << "Wrong rowDim in " << matrixNames[i]
			     << ", format " << ms.getFormat() << std::endl
			     << "Got " << m << ", should be " << rowDim << std::endl;
			fail = failThis = true;
		}
		if( !failThis && n != colDim ) {
			std::cout << "Wrong colDim in " << matrixNames[i]
			     << ", format " << ms.getFormat() << std::endl
			     << "Got " << n << ", should be " << colDim << std::endl;
			fail = failThis = true;
		}
		while( !failThis && ms.nextTriple(m,n,v) ) {
			if(!f.isZero(v)) --nzCount;
			if( m >= rowDim ) {
				std::cout << "Row index out of bounds in "
				     << matrixNames[i]
				     << ", format " << ms.getFormat() << std::endl
				     << "Got " << m << ", should be less than "
				     << rowDim << std::endl;
				fail = failThis = true;
				break;
			}
			if( n >= colDim ) {
				std::cout << "Column index out of bounds in "
				     << matrixNames[i]
				     << ", format " << ms.getFormat() << std::endl
				     << "Got " << n << ", should be less than "
				     << colDim << std::endl;
				fail = failThis = true;
				break;
			}
			if( matrix[m][n] != v ) {
				std::cout << "Invalid entry in "
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
			failThis = true;
		}
		if( !failThis && nzCount > 0 ) {
			std::cout << "Not enough entries in " << matrixNames[i]
			     << ", format " << ms.getFormat() << std::endl;
			fail = failThis = true;
		}
		if( !failThis && nzCount < 0 ) {
			std::cout << "Duplicate entries in " << matrixNames[i]
			     << ", format " << ms.getFormat() << std::endl;
			fail = failThis = true;
		}

		fin.seekg(0,std::ios::beg);
		MatrixStream<TestField> ms2(f,fin);
		std::vector<TestField::Element> array;
		ms2.getArray(array);
		if( array.size() != rowDim*colDim ) {
			std::cout << "Array given wrong size in " << matrixNames[i]
			     << ", format " << ms2.getFormat() << std::endl
			     << "Got " << array.size() << ", should be "
			     << rowDim*colDim << std::endl;
			fail = failThis = true;
		}
		for( m = 0; m < rowDim; ++m ) {
			for( n = 0; n < colDim; ++n ) {
				if( array[m*colDim+n] != matrix[m][n] ) {
					std::cout << "Invalid entry in getArray of "
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
			std::cout << "Test failed for " << matrixNames[i]
			     << ", format " << ms.getFormat() << std::endl
			     << "Error code: " << ms.getError()
			     << ", line number: " << ms.getLineNumber() << std::endl;
		}
	}

	if( 	testBlackBox< DenseMatrix<TestField> >
			( matrixNames[0], "Dense BlackBox Matrix" )
	  ) fail = true;
	if( 	testBlackBox< SparseMatrix<TestField> >
			( matrixNames[0], "Sparse BlackBox Matrix" )
	  ) fail = true;
	if( 	testBlackBox< BlasBlackbox<TestField> >
			( matrixNames[0], "BLAS BlackBox Matrix" )
	  ) fail = true;
	
	if( fail ) {
		std::cout << "FAIL: matrix-stream" << std::endl;
		return 1;
	}
	std::cout << "matrix-stream Passed" << std::endl;
	return 0;
}
