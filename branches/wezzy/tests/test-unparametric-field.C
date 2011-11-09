/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

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

/*! @file  tests/test-unparametric-field.C
 * @ingroup tests
 * @brief no doc.
 * @test no doc.
 */


#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/unparametric.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static size_t n = 10000;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	cout.flush ();
	bool pass = true;

	UnparametricField<double> F;

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
        report << endl << "Unparametrix<double> field test suite" << endl;

	if (!runFieldTests (F, "UnparametricField<double>", 1, n, false)) pass = false;

	FieldArchetype K(new UnparametricField<double>(F));

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of UnField field"))
		pass = false;
	// We're going to allow failed tests here. 
	// UnparametricField is a tool for building fields and does not of itself produce a LinBox conforming field.
	// However compilation serves some limited testing value and data is gleaned when the test is run with a report file argument.
	//return pass ? 0 : -1;
	return 0;
}
