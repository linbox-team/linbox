/* -*- mode: c; style: linux -*- */

/* tests/test-ntl-zz_p.cpp
 * Copyright (C) 2002 William J. Turner
 *
 * Written by William J. Turner <wjturner@math.ncsu.edu>
 *
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/ntl.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{

        static Argument args[] = {
                { '\0' }
        };

        parseArguments (argc, argv, args);

	cout << "UnparametricField<NTL::RR> field test suite" << endl << endl;
	cout.flush ();
	bool pass = true;

	UnparametricField<NTL::RR> F;

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	if (!testField<UnparametricField<NTL::RR> > (F, "Testing UnparametricField<NTL::RR> field"))
		pass = false;

#if 0
	FieldArchetype K(new UnparametricField<NTL::RR>);

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of UnField<NTL::RR> field"))
		pass = false;
#endif

	return pass ? 0 : -1;
}
