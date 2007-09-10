/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-hilbert.C  -bds
 */

#include "linbox/linbox-config.h"

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
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);
	typedef Modular<uint32> Field;
	Field F (q);

	srand (time (NULL));

	commentator.start("Hilbert matrix blackbox test suite", "Hilbert");

	typedef vector<Field::Element> Vector;
	typedef Hilbert<Field> BB;
	BB A (F, n);

	pass = pass && testBlackbox (F, A);

	commentator.stop("Hilbert matrix blackbox test suite");
	return pass ? 0 : -1;
}
