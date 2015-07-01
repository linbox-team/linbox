/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-quad-matrix.C
 * using generic testBlackbox  -bds
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/blackbox/quad-matrix.h"

#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	ofstream report;

	bool pass = true;

	static size_t n = 20;
	static integer q = 101;
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
	typedef ZOQuad <Field> BlackBox;

	Field F (q);
	Field::Element d; 
	F.init (d, 1);

	ScalarMatrix<Field> A (F, n, d); // a small identity.
        BlackBox AA(A);
	pass = pass && testBlackbox(F, AA);

  size_t *rows, *cols, i;
  rows = new size_t[3 * n + 1 - 2];
  cols = new size_t[3 * n + 1 - 2];

  // "arrow" matrix
  for(i = 0; i < n; i++) { rows[i] = 0; cols[i] = i; } // first row
  for(i = 0; i < n - 1; i++) { rows[n+2*i] = i + 1; cols[n+2*i] = 0; rows[n+2*i+1] = i + 1; cols[n+2*i+1] = i + 1; } // first col and the diag
	ZeroOne<Field> B(F, rows, cols, n, n, 3 * n - 2);
        BlackBox BB(B);
	pass = pass && testBlackbox(F, BB);

#if 1
	SideBySide<Field> C (AA, BB);
        BlackBox CC(C);
	pass = pass && testBlackbox(F, CC);

	OverUnder<Field> D(CC, CC);
        BlackBox DD(D); 
	pass = pass && testBlackbox(F, DD);
#endif

	return pass ? 0 : -1;
}

