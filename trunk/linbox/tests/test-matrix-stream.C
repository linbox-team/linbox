#include <iostream>
#include <fstream>
#include <string>
#include <linbox/field/unparametric.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/integer.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/sparse.h>

using std::cout;
using std::endl;
using namespace LinBox;

const int nMatrices = 7;
char* matrixNames[nMatrices] = {
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
	cout << "\tTesting " << BBName << endl;
	std::ifstream fin( filename );
	if( !fin ) {
		cout << "Could not open " << filename << endl;
		return true;
	}
	MatrixStream<TestField > ms(f, fin);
	BB m( ms );
	if( m.rowdim() != rowDim ) {
		cout << "Wrong rowDim in " << BBName << endl
		     << "Got " << m.rowdim() << ", should be " << rowDim << endl;
		return true;
	}
	if( m.coldim() != colDim ) {
		cout << "Wrong colDim in " << BBName << endl
		     << "Got " << m.coldim() << ", should be " << colDim << endl;
		return true;
	}
	bool fail = false;
	for( size_t i = 0; i < rowDim; ++i ) {
	    for( size_t j = 0; j < colDim; ++j ) {
	    	if( matrix[i][j] != m.getEntry( i, j ) ) {
			cout << "Invalid entry in " << BBName << " at index ("
			     << i << "," << j << ")" << endl
			     << "Got " << m.getEntry(i,j) << ", should be "
			     << matrix[i][j] << endl;
			fail = true;
		}
	    }
	}
	return fail;
}

int main() {
	bool fail = false;
	bool failThis;
	cout << "Testing matrix-stream..." << endl;

	for (int i = 0; i < 19; ++i) {matrix[2][2] *= 10; matrix[2][2] += 8; }

	for( int i = 0; i < nMatrices; ++i ) {
		failThis = false;
		std::ifstream fin(matrixNames[i]);
		if( !fin ) {
			cout << "Could not open " << matrixNames[i] << endl;
			fail = true;
			continue;
		}
		cout << "\tTesting " << matrixNames[i] << endl;
		MatrixStream<TestField > ms(f,fin);
		int nzCount = nonZeros;
		size_t m, n;
		integer v;
		if(!ms.getDimensions(m,n)) fail = failThis = true;
		if( !failThis && m != rowDim ) {
			cout << "Wrong rowDim in " << matrixNames[i]
			     << ", format " << ms.getFormat() << endl
			     << "Got " << m << ", should be " << rowDim << endl;
			fail = failThis = true;
		}
		if( !failThis && n != colDim ) {
			cout << "Wrong colDim in " << matrixNames[i]
			     << ", format " << ms.getFormat() << endl
			     << "Got " << n << ", should be " << colDim << endl;
			fail = failThis = true;
		}
		while( !failThis && ms.nextTriple(m,n,v) ) {
			if(!f.isZero(v)) --nzCount;
			if( m >= rowDim ) {
				cout << "Row index out of bounds in "
				     << matrixNames[i]
				     << ", format " << ms.getFormat() << endl
				     << "Got " << m << ", should be less than "
				     << rowDim << endl;
				fail = failThis = true;
				break;
			}
			if( n >= colDim ) {
				cout << "Column index out of bounds in "
				     << matrixNames[i]
				     << ", format " << ms.getFormat() << endl
				     << "Got " << n << ", should be less than "
				     << colDim << endl;
				fail = failThis = true;
				break;
			}
			if( matrix[m][n] != v ) {
				cout << "Invalid entry in "
				     << matrixNames[i]
				     << ", format " << ms.getFormat() << endl
				     << "Got " << v << " at index (" << m
				     << "," << n << "), should be "
				     << matrix[m][n] << endl;
				fail = failThis = true;
				break;
			}
		}
		if( ms.getError() != END_OF_MATRIX ) failThis = true;
		if( !failThis && nzCount > 0 ) {
			cout << "Not enough entries in " << matrixNames[i]
			     << ", format " << ms.getFormat() << endl;
			fail = failThis = true;
		}
		if( !failThis && nzCount < 0 ) {
			cout << "Duplicate entries in " << matrixNames[i]
			     << ", format " << ms.getFormat() << endl;
			fail = failThis = true;
		}
		if( failThis ) {
			cout << "Test failed for " << matrixNames[i]
			     << ", format " << ms.getFormat() << endl
			     << "Error code: " << ms.getError()
			     << ", line number: " << ms.getLineNumber() << endl;
		}
	}

	if( 	testBlackBox< DenseMatrix<TestField> >
			( matrixNames[0], "Dense BlackBox Matrix" )
	  ) fail = true;
	if( 	testBlackBox< SparseMatrix<TestField> >
			( matrixNames[0], "Sparse BlackBox Matrix" )
	  ) fail = true;
	
	if( fail ) {
		cout << "FAIL: matrix-stream" << endl;
		return 1;
	}
	cout << "matrix-stream Passed" << endl;
	return 0;
}
