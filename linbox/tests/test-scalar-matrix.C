/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-scalar-matrix.C
 * using generic testBlackbox  -bds
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/blackbox/scalar-matrix.h"

#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	ofstream report;

	bool pass = true;

	static size_t n = 20;
	static integer q = 2147483647U;
	static int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)",        TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 1)",          TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);

	srand (time (NULL));

	cout << endl << "Scalar matrix black box test suite" << endl;

	typedef Modular<uint32> Field;

	Field F (q);
	Field::Element d; 
	F.init (d, -1);

	typedef ScalarMatrix <Field> Blackbox;
	Blackbox A; // Test the default constructor
	pass = pass && testBlackbox(F, A);

	Blackbox B (F, n, d); // Test a small one.
	pass = pass && testBlackbox(F, B);

	Blackbox C (F, 100000, d); // Test a large one.
	pass = pass && testBlackbox(F, C);

	return pass ? 0 : -1;
}
