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
typedef GivaroPolyQuotient<PolyDom> QuotRing;
typedef GivaroPolyRandIter<PolyRing> PolyRandIter;
typedef typename PolyRing::Element PolyElement;
typedef MatrixDomain<PolyRing> MatDom;
typedef typename MatDom::OwnMatrix Mat;
typedef MatrixDomain<QuotRing> QuotMatDom;
typedef typename QuotMatDom::OwnMatrix QMat;
typedef IliopoulosDomain<PolyRing> IliopoulosDom;
typedef SmithFormLocal<QuotRing> LocalDom;
typedef DenseVector<PolyRing> FactorVector;

void local(Mat &M, PolyElement &f, QuotRing &R) {
	size_t n = M.coldim();
	QMat A(R, n, n);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			PolyElement tmp;
			M.getEntry(tmp, i, j);
			A.setEntry(i, j, tmp);
		}
	}
	
	cout << "Starting Local" << endl;
	LocalDom sfl;
	list<PolyElement> l;
	sfl(l, A, R);
	
	list<PolyElement>::const_iterator iterator;
	for (iterator = l.begin(); iterator != l.end(); ++iterator) {
		R.write(cout, *iterator) << endl;
	}
	cout << endl;
}

void ilio(Mat &M, PolyElement &f, MatDom &MD, PolyRing &R) {
	IliopoulosDom sfi(R);
	
	size_t n = M.coldim();
	Mat A(MD.field(), n, n);
	MD.copy(A, M);
	
	FactorVector l;
	l.resize(n);
	sfi.smithForm(l, A, f);
	
	for (size_t i = 0; i < l.size(); i++) {
		R.write(cout, l[i]) << endl;
	}
	cout << endl;
}

void generateM(Mat &M, MatDom &MD, PolyRing &R, PolyElement &f, PolyElement &g) {
	size_t n = 3;
	
	integer max;
	R.convert(max, f);
	max *= 10;
	
	Mat L(R, n, n);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < i; j++) {
			PolyElement e;
			R.init(e, rand() % max);
			R.modin(e, f);
			L.setEntry(i, j, e);
		}
	}
	for (size_t i = 0; i < n; i++) {
		L.setEntry(i, i, R.one);
	}
	
	Mat T(R, n, n);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = i+1; j < n; j++) {
			PolyElement e;
			R.init(e, rand() % max);
			R.modin(e, f);
			T.setEntry(i, j, e);
		}
	}
	for (size_t i = 0; i < n; i++) {
		T.setEntry(i, i, R.one);
	}
	
	Mat D(R, n, n);
	D.setEntry(0, 0, R.one);
	D.setEntry(1, 1, g);
	PolyElement tmp;
	R.mul(tmp, g, g);
	R.mulin(tmp, g);
	D.setEntry(2, 2, tmp);
	
	MD.mul(M, L, D);
	MD.mulin(M, T);
	
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			PolyElement e;
			M.getEntry(e, i, j);
			R.write(cout << i << ", " << j << ": ", e) << endl;
		}
	}
}

int main(int argc, char* argv[]) {
	size_t p = 3;
	size_t n = 4;
	
	Field F(p);
	PolyDom PD(F, "x");
	PolyRing R(PD);
	PolyRandIter PRI(R, 0, 10);
	MatDom MD(R);
	
	PolyElement g, f;
	R.init(g, 13);
	R.write(cout, g) << endl;
	R.mul(f, g, g);
	R.mulin(f, f);
	R.write(cout, f) << endl;
	
	QuotRing QR(PD, f);
	
	Mat M(R, n, n);
	generateM(M, MD, R, f, g);
	
	cout << endl;
	
	ilio(M, f, MD, R);
	
	local(M, f, QR);
}