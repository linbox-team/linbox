#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
//#include <omp.h>

//#define LINBOX_USES_OMP 1
#include "linbox/ring/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include <linbox/ring/polynomial-ring.h>
#include "linbox/algorithms/smith-form-kannan-bachem.h"

using namespace std;
using namespace LinBox;

int main(int argc, char **argv) {
	typedef Givaro::Modular<double> Field;
	typedef PolynomialRing<Field> PolyDom;
	typedef PolynomialRing<PolyDom> PolyRing;
	typedef PolyRing::Element Element;
	typedef DenseVector<PolyRing> PolyVector;
	typedef MatrixDomain<PolyRing> PolyMatDom;
	typedef typename PolyMatDom::OwnMatrix Matrix;
	
	int p = 3;
	int n = 3;
	
	Field F(p);
	PolyDom PD(F,"x");
	PolyRing R(PD);
	PolyMatDom PMD(R);
	
	Matrix M(R, n, n);
	
	Element a; // x + 2
	R.init(a, 5);
	
	Element b; // 2
	R.init(b, 2);
	
	Element c; // x
	R.init(c, 3);
	
	M.setEntry(0,0,a);
	M.setEntry(0,1,b);
	M.setEntry(0,2,b);
	
	M.setEntry(1,1,c);
	
	M.setEntry(2,2,c);
	
	SmithFormKannanBachemDomain<PolyMatDom> SFKB(PMD);
	
	PolyVector factors(R);
	factors.resize(n);
	
	SFKB.solve(factors, M);
	
	R.init(a, 1);
	R.init(b, 3);
	R.init(c, 15);
	
	bool pass = true;
	
	pass = pass && R.areEqual(a, factors[0]);
	pass = pass && R.areEqual(b, factors[1]);
	pass = pass && R.areEqual(c, factors[2]);
	
	return pass ? 0 : -1;
}
