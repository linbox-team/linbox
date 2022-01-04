#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
//#include <omp.h>
//#define LINBOX_USES_OMP 1

#include "linbox/ring/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "givaro/givpoly1.h"
#include "linbox/ring/polynomial-ring.h"
#include "linbox/algorithms/smith-form-kannan-bachem.h"

using namespace std;
using namespace LinBox;

int main(int argc, char **argv) {
	typedef Givaro::Modular<double> Field;
	typedef PolynomialRing<Field> PolyDom;
	typedef PolyDom::Element Element;
	typedef std::vector<Element> PolyVector;
	typedef MatrixDomain<PolyDom> PolyMatDom;
	typedef typename PolyMatDom::OwnMatrix Matrix;

	int p = 3;
	int n = 3;

	Field F(p);
	PolyDom R(F,'x');
	PolyMatDom PMD(R);

	Matrix M(R, n, n);

	Element a; // x + 2
    R.init(a, {2,1} );

	Element b; // 2
	R.init(b, 2);

	Element c; // x
	R.init(c, {0,1} );

	M.setEntry(0,0,a);
	M.setEntry(0,1,b);
	M.setEntry(0,2,b);

	M.setEntry(1,1,c);

	M.setEntry(2,2,c);

    M.write(std::clog << "smith(", Tag::FileFormat::linalg) << ",x)";

	SmithFormKannanBachemDomain<PolyDom> SFKB(R);

	PolyVector factors; factors.reserve(n);

	SFKB.solve(factors, M);

    std::clog << " is |";
    for(const auto& iter: factors) R.write(std::clog,iter) << '|';


	bool pass = true;

	R.init(a, 1);		// 1
	R.init(b, {0,1});	// x
	R.init(c, {0,2,1});	// 2x+x^2
	pass = pass
        && R.areEqual(a, factors[0])
        && R.areEqual(b, factors[1])
        && R.areEqual(c, factors[2]);

    if (!pass) {
        std::cerr << " *** ERROR ***" << std::endl;
        R.write(std::cerr << "1: ", a) << std::endl;
        R.write(std::cerr << "x: ", b) << std::endl;
        R.write(std::cerr << "0+x(2+x): ", c) << std::endl;
        R.write(std::cerr << "f[0]: ", factors[0]) << std::endl;
        R.write(std::cerr << "f[1]: ", factors[1]) << std::endl;
        R.write(std::cerr << "f[2]: ", factors[2]) << std::endl;
    } else std::clog << ", PASSED." << std::endl;
        

	return pass ? 0 : -1;
}
