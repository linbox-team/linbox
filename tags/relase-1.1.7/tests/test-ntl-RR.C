/* tests/test-ntl-zz_p.cpp
 * Copyright (C) 2002 William J. Turner
 * Copyright (C) LinBox
 *
 * Written by William J. Turner <wjturner@math.ncsu.edu>
 *
 * see COPYING file for license
 *
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/ntl.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static size_t n = 10000;
	static int iterations = 1;

        static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
                { '\0' }
        };

        parseArguments (argc, argv, args);

	commentator.start("UnparametricField<NTL::RR> field test suite", "UnparametricField<NTL::RR>");
	bool pass = true;

	UnparametricField<NTL::RR> F;

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	if (!runFieldTests (F, "UnparametricField<NTL::RR>", iterations, n, false)) pass = false;

#if 0
	FieldArchetype K(new UnparametricField<NTL::RR>);

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of UnField<NTL::RR> field"))
		pass = false;
#endif

	commentator.stop("UnparametricField<NTL::RR> field test suite");
	// We're going to allow failed tests here because the floating-point
	// approximation tends to screw things up anyway.

	// -bds:  Well, compilation is checked at least.  Work needed: A meaningful test is falsifyable.

	return 0;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
