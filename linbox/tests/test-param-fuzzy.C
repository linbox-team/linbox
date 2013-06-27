
/* tests/test-param-fuzzy.C
 * Copyright (C) 2002 David Saunders
 * shamelessly mutated from one of the other field tests.
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


/*! @file  tests/test-param-fuzzy.C
 * @ingroup tests
 * @brief  no doc
 * @test no doc
 */



#include "linbox-config.h"

#include <iostream>
#include <fstream>


#include "linbox/field/param-fuzzy.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	//static integer q = 2147483647U;
	static size_t n = 10000;
	static unsigned int iterations = 10;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	bool pass = true;

	ParamFuzzy F(.0000001);

	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
        report << endl << "ParamFuzzy field approximation test suite" << endl;

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

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

