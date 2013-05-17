
/* tests/test-ntl-lzz_p.cpp
 * Copyright (C) 2002 William J. Turner
 *
 * Written by William J. Turner <wjturner@math.ncsu.edu>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file  tests/test-ntl-lzz_p.C
 * @ingroup tests
 * @brief no doc.
 * @test no doc.
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
        static integer q = 1073741789;
	static size_t n = 10000;
	static int iterations = 1;

        static Argument args[] = {
                { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
        };

        parseArguments (argc, argv, args);

	commentator().start("NTL_zz_p field test suite", "NTL_zz_p");
	bool pass = true;

	//NTL::zz_p::init(q);
	//UnparametricField<NTL::zz_p> F(q);
	// NTL_zz_p F(q);
	NTL_zz_p F;

	if (F.characteristic() != q)
		return false ;

	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	if (!runFieldTests (F, "NTL_zz_p", (unsigned int)iterations, n, false)) pass = false;

#if 0
	FieldArchetype K(new NTL_zz_p(101));

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of UnField<NTL::zz_p> field"))
		pass = false;
#endif

	commentator().stop("NTL_zz_p field test suite");
	return pass ? 0 : -1;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

