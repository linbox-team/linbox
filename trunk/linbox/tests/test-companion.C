/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-scalar-matrix.C
 * using generic testBlackbox  -bds
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/blackbox/companion.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	ofstream report;

	bool pass = true;

	static size_t n = 10;
	static integer q = 2147483647U;
	static int iterations = 10;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)",        TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",          TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);

	srand (time (NULL));

	cout << endl << "Companion matrix black box test suite" << endl;

	typedef Modular<uint32> Field;
	typedef vector <Field::Element> Vector;
	typedef vector <Field::Element> Polynomial;
	typedef Companion<Field> Blackbox;

	Field F (q);
	Field::Element d; 
	F.init (d, -1);
	Polynomial p(n+1, d);  
	F.init (d, 1); F.assign(p[n], d);

	Blackbox A (F, p);

	pass = pass && testBlackbox<Field, Blackbox>(F, A);

	Blackbox B (F, n);

	pass = pass && testBlackbox<Field, Blackbox>(F, B);

	return pass ? 0 : -1;
}
