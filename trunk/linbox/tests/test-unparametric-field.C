/* -*- mode: c; style: linux -*- */

/* tests/test-unparametric-field.C
 * Copyright (C) 2002 William J. Turner
 *
 * Written by William J. Turner <wjturner@math.ncsu.edu>
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Renamed from test-unparametric-field.cpp to test-unparametric-field.C, so
 * that we are using the same file naming conventions thoughout the library.
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/unparametric.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static integer q = 2147483647U;

	static Argument args[] = {
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	cout << "Unparametrix<double> field test suite" << endl << endl;
	cout.flush ();
	bool pass = true;

	UnparametricField<double> F;

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	if (!testField<UnparametricField<double> > (F, "Testing UnparametricField<double> field"))
		pass = false;

#if 0
	FieldArchetype K(new UnparametricField(101));

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of UnField field"))
		pass = false;
#endif

	return pass ? 0 : -1;
}
