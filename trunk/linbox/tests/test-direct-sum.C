/* tests/test-direct-sum.C
 * Copyright (C) LinBox
 * Written by David Saunders
 *
 * See COPYING for license information.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/direct-sum.h"

#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 101;
	static int iterations1 = 100;
	static int iterations2 = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.",   TYPE_INT,     &iterations1 },
		{ 'j', "-j J", "Apply test matrix to J vectors.",         TYPE_INT,     &iterations2 },
		{ '\0' }
	};

	typedef Modular<uint32> Field;
	typedef vector<Field::Element> Vector;

	parseArguments (argc, argv, args);
	Field F (q);
	Field::Element k;

	commentator.start("DirectSum black box test suite", "direct sum");
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	F.init(k, 5);
	ScalarMatrix<Field> B(F, 10, k);
	F.init(k, 2);
	ScalarMatrix<Field> C(F, 5, k);
	
	DirectSum<ScalarMatrix<Field>, ScalarMatrix<Field> > A(&B, &C);
	pass = pass && testBlackbox(A);
	DirectSum<ScalarMatrix<Field>, ScalarMatrix<Field> > D(B, C);
	pass = pass && testBlackbox(D);

	commentator.stop("DirectSum black box test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
