
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
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file  tests/test-unparametric-field.C
 * @ingroup tests
 * @brief no doc.
 * @test no doc.
 */


#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>


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

	//cout.flush ();
	bool pass = true;

	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
        report << endl << "Unparametrix<double> field test suite" << endl;

	UnparametricField<double> F;

	if (!runFieldTests (F, "UnparametricField<double>", 1, n, false)) pass = false;

// archetype no longer works
#if 0
	UnparametricField<double> * F1 = new UnparametricField<double>(F);
	FieldArchetype K(F1);
	delete F1;

	// FieldArchetype K(new UnparametricField<double>(F);) // leaks !!!

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of UnField field"))
		pass = false;
#endif
	// We're going to allow failed tests here.
	// UnparametricField is a tool for building fields and does not of itself produce a LinBox conforming field.
	// However compilation serves some limited testing value and data is gleaned when the test is run with a report file argument.
	// return 0;
	// but for now accept the test.
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
