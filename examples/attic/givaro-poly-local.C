#include <linbox/linbox-config.h>

#include <iostream>
#include <string>
#include <vector>
#include <list>

#include <linbox/util/timer.h>

#include <linbox/ring/modular.h>

#include <linbox/ring/pir-modular-int32.h>
#include <linbox/ring/givaro-poly.h>
#include <linbox/ring/givaro-poly-local.h>
#include <linbox/ring/givaro-poly-quotient.h>
#include <linbox/randiter/givaro-poly.h>

using namespace LinBox;
using namespace std;

typedef Givaro::Modular<double> Field;
typedef typename Field::Element Element;
typedef Givaro::Poly1Dom<Field,Givaro::Dense> PolyDom;
typedef GivaroPoly<Field> PolyRing;

typedef GivaroPolyQuotient<PolyDom> QuotRing;
typedef typename QuotRing::Element QElm;

typedef GivaroPolyLocal<PolyDom> LocalRing;
typedef typename LocalRing::Element LElm;

typedef GivaroPolyRandIter<PolyRing> PolyRandIter;
typedef typename PolyRing::Element PolyElement;

bool testMul(QuotRing QR, LocalRing LR, QElm qa, QElm qb, LElm la, LElm lb) {
	QElm qc;
	QR.mul(qc, qa, qb);

	LElm lc;
	LR.mul(lc, la, lb);

	if (LR.areEqual(lc, qc)) {
		return true;
	}

	QR.write(cout << "(", qa) << ") * (";
	QR.write(cout, qb) << ") = ";
	QR.write(cout, qc) << endl;

	LR.write(cout << "(", la) << ") * (";
	LR.write(cout, lb) << ") = ";
	LR.write(cout, lc) << endl << endl;

	LR.write2(cout, la) << endl;
	LR.write2(cout, lb) << endl;

	return false;
}

bool testAdd(QuotRing QR, LocalRing LR, QElm qa, QElm qb, LElm la, LElm lb) {
	QElm qc;
	QR.add(qc, qa, qb);

	LElm lc;
	LR.add(lc, la, lb);

	if (LR.areEqual(lc, qc)) {
		return true;
	}

	QR.write(cout << "(", qa) << ") + (";
	QR.write(cout, qb) << ") = ";
	QR.write(cout, qc) << endl;

	LR.write(cout << "(", la) << ") + (";
	LR.write(cout, lb) << ") = ";
	LR.write(cout, lc) << endl << endl;

	LR.write2(cout, la) << endl;
	LR.write2(cout, lb) << endl;

	return false;
}

bool testAxpy(
	QuotRing QR, LocalRing LR,
	QElm qa, QElm qb, QElm qc,
	LElm la, LElm lb, LElm lc
) {
	QElm qd;
	QR.axpy(qd, qa, qb, qc);

	LElm ld;
	LR.axpy(ld, la, lb, lc);

	if (LR.areEqual(ld, qd)) {
		return true;
	}

	QR.write(cout << "(", qa) << ") * (";
	QR.write(cout, qb) << ") + (";
	QR.write(cout, qc) << ") = ";
	QR.write(cout, qd) << endl;

	LR.write(cout << "(", la) << ") * (";
	LR.write(cout, lb) << ") + (";
	LR.write(cout, lc) << ") = ";
	LR.write(cout, ld) << endl;

	LR.write2(cout, la) << endl;
	LR.write2(cout, lb) << endl;
	LR.write2(cout, lc) << endl;

	return false;
}

int main(int argc, char* argv[]) {
	size_t p = 3;
	size_t e = 3;

	Field F(p);
	PolyRing R(F, "x");
	PolyRandIter PRI(R, 0, 10);

	PolyElement g, f;
	R.init(g, {2, 1, 1});
	R.write(cout << "irred: ", g) << endl;

	R.pow(f, g, e);
	R.write(cout << "quotient: ", f) << endl;

	integer max;
	R.convert(max, f);

	QuotRing QR(R, f);
	LocalRing LR(R, g, e);

	QElm qa, qb, qc;
	LElm la, lb, lc;

	for (integer i = 1; i <= max; i++) {
		QR.init(qa, i);
		LR.init(la, i);

		for (integer j = 1; j <= max; j++) {
			cout << i << ", " << j << ":" << endl;

			QR.init(qb, j);
			LR.init(lb, j);

			//if (!QR.isDivisor(qa, qb)) {
			//	continue;
			//}

			if (!testMul(QR, LR, qa, qb, la, lb)) {
				return 0;
			}

			if (!testAdd(QR, LR, qa, qb, la, lb)) {
				return 0;
			}

			for (integer k = 1; k <= max; k++) {
				QR.init(qc, k);
				LR.init(lc, k);

				if (!testAxpy(QR, LR, qa, qb, qc, la, lb, lc)) {
					return 0;
				}
			}
		}
	}
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
