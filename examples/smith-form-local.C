#include <linbox/linbox-config.h>

#include <iostream>
#include <string>
#include <vector>
#include <list>

using namespace std;


#include <linbox/ring/modular.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/algorithms/smith-form-sparseelim-local.h>

#include <linbox/util/timer.h>

#include <linbox/ring/local2_32.h>
#include <linbox/ring/pir-modular-int32.h>
#include <linbox/ring/givaro-poly.h>
#include <linbox/ring/givaro-poly-quotient.h>
#include <linbox/randiter/givaro-poly.h>
#include <linbox/algorithms/smith-form-local.h>
#include <linbox/algorithms/smith-form-iliopoulos2.h>
#include <linbox/algorithms/smith-form-adaptive.h>

using namespace LinBox;
using namespace std;

typedef Givaro::Modular<double> Field;
typedef typename Field::Element Element;
typedef Givaro::Poly1Dom<Field,Givaro::Dense> PolyDom;
typedef GivaroPoly<PolyDom> PolyRing;
typedef GivaroPolyRandIter<PolyRing> PolyRandIter;
typedef typename PolyRing::Element PolyElement;
typedef MatrixDomain<PolyRing> MatDom;
typedef typename MatDom::OwnMatrix Mat;
typedef IliopoulosDomain<PolyRing> IliopoulosDom;
typedef DenseVector<PolyRing> FactorVector;

void local(list<PolyElement> &l, Mat &M, PolyElement &f, MatDom &MD, PolyDom &R) {
	typedef MatrixDomain<GivaroPolyQuotient<PolyDom>> QuotMatDom;
	typedef typename QuotMatDom::OwnMatrix QMat;
	
	GivaroPolyQuotient<PolyDom> QR(R, f);
	QuotMatDom QMD(QR);
	
	size_t n = M.coldim();
	QMat A(QR, n, n);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			PolyElement tmp;
			M.getEntry(tmp, i, j);
			A.setEntry(i, j, tmp);
		}
	}
	
	SmithFormLocal<GivaroPolyQuotient<PolyDom>> sfl;
	sfl(l, A, QR);
}

void ilio(FactorVector &l, Mat &M, PolyElement &f, MatDom &MD, PolyRing &R) {
	IliopoulosDom sfi(R);
	
	size_t n = M.coldim();
	Mat A(MD.field(), n, n);
	MD.copy(A, M);
	sfi.smithForm(l, A, f);
}

int main(int argc, char* argv[]) {
	size_t p = 3;
	size_t n = 3;
	
	Field F(p);
	PolyDom PD(F, "x");
	PolyRing R(PD);
	PolyRandIter PRI(R, 0, 10);
	MatDom MD(R);
	
	PolyElement f;
	R.init(f, 13);
	R.write(cout, f) << endl;
	R.mulin(f, f);
	R.mulin(f, f);
	R.write(cout, f) << endl;
	
	integer max;
	R.convert(max, f);
	max *= 10;
	
	Mat M(R, n, n);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			PolyElement e;
			R.init(e, rand() % max);
			R.modin(e, f);
			R.write(cout, e) << endl;
			M.setEntry(i, j, e);
		}
	}
	
	list<PolyElement> l;
	local(l, M, f, MD, PD);
	
	cout << endl;
	list<PolyElement>::const_iterator iterator;
	for (iterator = l.begin(); iterator != l.end(); ++iterator) {
		R.write(cout, *iterator) << endl;
	}
	
	FactorVector fl;
	fl.resize(n);
	ilio(fl, M, f, MD, R);
	
	cout << endl;
	for (size_t i = 0; i < fl.size(); i++) {
		R.write(cout, fl[i]) << std::endl;
	}
}