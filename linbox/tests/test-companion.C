/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-scalar-matrix.C
 * using generic testBlackbox  -bds
 */

#include "linbox/linbox-config.h"

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
	bool pass = true;

	static size_t n = 10;
	static integer q = 2147483647U;
	static int iterations = 1; // was 10

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.",        TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.",          TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	srand (time (NULL));

	commentator.start("Companion matrix black box test suite", "companion");

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

	pass = pass && testBlackbox(A);

	Blackbox B (F, n);

	pass = pass && testBlackbox(B);

	commentator.stop("companion matrix black box test suite");

	return pass ? 0 : -1;
}
