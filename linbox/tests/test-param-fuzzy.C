/* -*- mode: C++; style: linux -*- */

/* tests/test-param-fuzzy.C
 * Copyright (C) 2002 David Saunders
 * shamelessly mutated from one of the other field tests.
 *
 * See COPYING for license information.
 */

#include "linbox-config.h"

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

	static Argument args[] = {
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	cout << "ParamFuzzy field approximation test suite" << endl << endl;
	cout.flush ();
	bool pass = true;

	ParamFuzzy F(.0000001);

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	/* I am distressed that this field passes the testField()
	   We need a test that distinguishes exact fields from 
	   approximate ones.  -bds */

	if (!testField<ParamFuzzy > (F, "Testing ParamFuzzy field"))
		pass = false;

#if 0
	FieldArchetype K(new ParamFuzzy(.000001));

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of UnField field"))
		pass = false;
#endif

	DoubleRealApproximation G(.00000000001);
	if (!testField<DoubleRealApproximation > (F, "Testing DoubleRealApproximation field"))
		pass = false;

	return pass ? 0 : -1;

}
