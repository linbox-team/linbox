#include "linbox/linbox-config.h"

#include <iostream>

#include "givaro/givtimer.h"

#include "linbox/ring/modular.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/random-matrix.h"

using namespace LinBox;

typedef Givaro::GFqDom<int64_t> Field;
//typedef Givaro::Modular<double> Field;
typedef typename Field::RandIter RandIter;
typedef BlasMatrixDomain<Field> MatrixDom;
typedef BlasMatrix<Field> Matrix;

int main() {
	size_t p = 3, e = 13;
	
	Field F(p, e);
	//Field F(p);
	
	MatrixDom MD(F);
	
	Matrix A(F, 5, 5);
	Matrix B(F, 5, 5);
	Matrix C(F, 5, 5);
	
	RandIter RI(F);
	RandomDenseMatrix<RandIter, Field> RDM(F, RI);
	
	RDM.random(A);
	RDM.random(B);
	
	MD.mul(C, A, B);
	
	A.write(std::cout) << std::endl;
	B.write(std::cout) << std::endl;
	std::cout << "==" << std::endl;
	C.write(std::cout) << std::endl;
	
	return 0;
}