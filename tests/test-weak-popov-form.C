#include <linbox/linbox-config.h>
#include <vector>
#include <utility>

#include "linbox/ring/ntl.h"
#include "linbox/algorithms/weak-popov-form.h"

//#define LINBOX_USES_OMP 1
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/random-matrix.h"

using namespace LinBox;

typedef NTL_zz_pX PolynomialRing;
typedef typename PolynomialRing::Element Polynomial;

typedef MatrixDomain<PolynomialRing> PolyMatrixDom;
typedef typename PolyMatrixDom::OwnMatrix Matrix;
typedef WeakPopovFormDomain<PolynomialRing> WeakPopovFormDom;

bool writing = false;

int main(int argc, char* argv[]) {

	// text is written to cout iff a command line arg is present	
	if (argc > 1) writing = true;

	uint64_t p = 3;
	PolynomialRing R(p);
	PolyMatrixDom MD(R);
	WeakPopovFormDom PFD(R);
	
	Matrix M(R, 4, 4);
	M.setEntry(0, 1, R.one);
	M.setEntry(1, 1, R.one);
	M.setEntry(2, 1, R.one);
	M.setEntry(2, 2, R.one);
	
	Matrix V(R, 4, 1);
	V.setEntry(0, 0, R.one);
	V.setEntry(1, 0, R.one);
	V.setEntry(2, 0, R.one);
	V.setEntry(3, 0, R.one);
	
	if (writing) PFD.printMatrix(M);
	if (writing) std::cout << std::endl;
	
	std::vector<long> pivots;
	PFD.findPivots(pivots, M);
	
	if (writing) std::cout << pivots << std::endl;
	
	std::pair<size_t, size_t> match;
	bool matchFound = PFD.findMatchingPivots(match, pivots);
	
	if (writing) std::cout << (matchFound ? "true: " : "false: ") << match << std::endl;
	
	PFD.extendedSolve(M, V);
	
	if (writing) PFD.printMatrix(M);
	if (writing) std::cout << std::endl;
	
	if (writing) PFD.printMatrix(V);
	if (writing) std::cout << std::endl;
	
	Polynomial det;
	PFD.solveDet(det, M);
	
	return 0;
}
