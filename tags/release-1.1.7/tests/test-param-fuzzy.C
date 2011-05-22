
/* tests/test-param-fuzzy.C
 * Copyright (C) 2002 David Saunders
 * shamelessly mutated from one of the other field tests.
 *
 * See COPYING for license information.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/param-fuzzy.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	//static integer q = 2147483647U;
	static size_t n = 10000;
	static int iterations = 10;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	cout << endl << "ParamFuzzy field approximation test suite" << endl;
	cout.flush ();
	bool pass = true;

	ParamFuzzy F(.0000001);

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	/* I am distressed that this field passes the testField()
	   We need a test that distinguishes exact fields from 
	   approximate ones.  -bds */

	if (!runFieldTests (F, "ParamFuzzy", iterations, n, false)) pass = false;

#if 0
	FieldArchetype K(new ParamFuzzy(.000001));

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of UnField field"))
		pass = false;
#endif

	if (!testField<DoubleRealApproximation > (F, "Testing DoubleRealApproximation field"))
		pass = false;

	return pass ? 0 : -1;

}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
