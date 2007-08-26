/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-scalar-matrix.C
 * using generic testBlackbox  -bds
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/blackbox/triplesbb.h"
#include "linbox/blackbox/transpose.h"

#include "test-common.h"
#include "test-generic.h"


template<class Vec>
bool eqVec(const Vec& a, const Vec& b)
{       bool good = true;
	for (typename Vec::size_type i = 0; i < a.size(); ++i) good = good && (a[i] == b[i]);
	return good;
}
using namespace LinBox;

int main (int argc, char **argv)
{
	ofstream report;

	bool pass = true;

	static size_t n = 200;
	static integer q = 2147483647U;
	static int iterations = 10;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	srand (time (NULL));

	commentator.start("triplesbb black box test suite", "triplesbb");

	typedef Modular<uint32> Field;
	typedef Field::Element Element;
	typedef vector <Element> Vector;
	typedef TriplesBB<Field> Blackbox;

	Field F (q);
	Element d; 
	F.init (d, -1);

	// set up the matrix
	std::vector<Element> values; 
	std::vector<size_t> rowP; 
	std::vector<size_t> colP; 
	for(int i = 1; i < 12; ++i)
	{
		values.push_back(i);
		rowP.push_back(i%n + 1);
		colP.push_back(i/n + 1);
	}

	Blackbox A(F, values, rowP, colP, n, n);

	pass = pass && testBlackbox<Field,Blackbox>(F, A);

	Transpose<Blackbox> B(&A);

	pass = pass && testBlackbox<Field,Transpose<Blackbox> >(F, B);

	Vector x(n), y(n), z(n);
	for(int i = 0; i < 5; ++i) x[i] = i;

	A.apply(y, x);
	B.applyTranspose(z, x);
	pass = pass && eqVec(y, z);

	B.apply(y, x);
	A.applyTranspose(z, x);
	pass = pass && eqVec(y, z);

	Blackbox C(F, n, n);
	for(size_t i = 0; i < rowP.size(); ++i) C.addEntry(values[i], rowP[i], colP[i]);
	pass = pass && testBlackbox<Field, Blackbox>(F, C);

	commentator.stop("triplesbb black box test suite");
	return pass ? 0 : -1;
}
