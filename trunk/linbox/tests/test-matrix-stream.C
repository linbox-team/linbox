#include <iostream>
#include <fstream>
#include <string>
#include <linbox/field/unparametric.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/integer.h>

using std::cout;
using std::endl;
using namespace LinBox;

const int nMatrices = 7;
char* matrixNames[nMatrices] = {"data/matrix-market-array.matrix",
				"data/maple-sparse1.matrix",
                                "data/maple-dense1.matrix",
                                "data/generic-dense.matrix",
				"data/sparse-row.matrix",
				"data/sms.matrix",
				"data/matrix-market-coordinate.matrix"};

const int rowDim = 11;
const int colDim = 11;
int nonZeros = 33;
integer matrix[rowDim][colDim] = {
                           {0, 0, 2, 3, 0, 0, 0, 0, 0, 1, 0},
                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			   {2, 0, 8888888888888888888, 1, -1, 0, 0, 0, 0, 0, 6},
			   {3, 0, 1, 4, 0, 0, 12, 0, 0, -13, 0},
			   {0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0},
			   {0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0},
			   {0, 0, 0, 12, 0, 0, 0, 0, 0, 0, 0},
			   {0, 0, 0, 0, 0, 1, 0, 500, 400, 300, 200},
			   {0, 0, 0, 0, 0, 0, 0, 400, 0, 0, 0},
			   {1, 0, 0, -13, 0, 1, 0, 300, 0, 10, 1},
			   {0, 0, 6, 0, 0, 0, 0, 200, 0, 1, 0} };

int main() {
	bool fail = false;
	bool failThis;
	cout << "Testing matrix-stream..." << endl;
	UnparametricField<integer> f;
	for( int i = 0; i < nMatrices; ++i ) {
		failThis = false;
		std::ifstream fin(matrixNames[i]);
		if( !fin ) {
			cout << "Could not open " << matrixNames[i] << endl;
			fail = true;
			continue;
		}
		MatrixStream<UnparametricField<integer> > ms(f,fin);
		int nzCount = nonZeros;
		int m, n;
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
		if( !failThis && nzCount > 0 ) {
			cout << "Not enough entries in " << matrixNames[i]
			     << ", format " << ms.getFormat() << endl;
			fail = failThis = true;
		}
		if( !failThis && nzCount < 0 ) {
			cout << "Warning: nonZeros value is incorrect.  "
			     << "test-matrix-stream.C needs fixing." << endl;
		}
		if( failThis ) {
			cout << "Test failed for " << matrixNames[i]
			     << ", format " << ms.getFormat() << endl
			     << "Error code: " << ms.getError()
			     << ", line number: " << ms.getLineNumber() << endl;
		}
		ms.~MatrixStream();
	}
	if( fail ) return 1;
	cout << "matrix-stream Passed" << endl;
	return 0;
}
