/* -*- mode: c; style: linux -*- */

/* tests/test-ntl-zz_p.cpp
 * Copyright (C) 2002 William J. Turner
 *
 * Written by William J. Turner <wjturner@math.ncsu.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/ntl-zz_p.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
        static integer q = 1073741789;
	static size_t n = 10000;
	static int iterations = 10;

        static Argument args[] = {
                { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 1073741789)", TYPE_INTEGER, &q },
		{ 'n', "-n N", "Set dimension of test vectors to NxN (default 10000)",      TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",           TYPE_INT,     &iterations },
                { '\0' }
        };

        parseArguments (argc, argv, args);

	cout << "UnparametricField<NTL::zz_p> field test suite" << endl << endl;
	cout.flush ();
	bool pass = true;

	//NTL::zz_p::init(q);
	//UnparametricField<NTL::zz_p> F(q);
	NTL_zz_p_Field F(q); 

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	if (!runFieldTests (F, "UnparametricField<NTL::zz_p>", iterations, n, false)) pass = false;

#if 0
	FieldArchetype K(new UnparametricField<NTL::zz_p>(101));

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of UnField<NTL::zz_p> field"))
		pass = false;
#endif

	return pass ? 0 : -1;
}
