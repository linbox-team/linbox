/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-hilbert.C  -bds
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/blackbox/hilbert.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

/* run generic testBlackbox on a Hilbert matrix */
int main (int argc, char **argv)
{
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
	typedef Modular<uint32> Field;
	Field F (q);

	srand (time (NULL));

	cout << endl << "Hilbert matrix blackbox test suite" << endl;

	typedef vector<Field::Element> Vector;
	typedef Hilbert<Field> BB;
	BB A (F, n);

	pass = pass && testBlackbox (F, A);

	return pass ? 0 : -1;
}
