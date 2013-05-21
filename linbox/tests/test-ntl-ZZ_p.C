
/* tests/test-ntl-ZZ_p.C
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


/*! @file  tests/test-ntl-ZZ_p.C
 * @ingroup tests
 * @brief  no doc
 * @test no doc
 */



#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>


#include "linbox/field/ntl.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
        static long q = 1073741789;
	static size_t n = 10000;
	static int iterations = 1;

        static Argument args[] = {
                { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
        };

        parseArguments (argc, argv, args);

	commentator().start("NTL_ZZ_p field test suite", "NTL_ZZ_p");
	bool pass = true;

	NTL::ZZ_p::init(NTL::to_ZZ(q));
	// UnparametricField<NTL::ZZ_p> F;
	// NTL_ZZ_p F(q); // XXX is there a q ?
	NTL_ZZ_p F;

	if (F.characteristic() != q)
		return false ;


	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	if (!runFieldTests (F, "NTL_ZZ_p", (unsigned int)iterations, n, false)) pass = false;

#if 0
	FieldArchetype K(new NTL_ZZ_p(101));

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of UnField<NTL::ZZ_p> field"))
		pass = false;
#endif

	// Testing big ints
	//
        commentator().start ("\t--Testing init/convert match");
        bool part_pass = true;
	// UnparametricField<NTL::ZZ_p> G;
	// NTL_ZZ_p G(q);
	NTL_ZZ_p G;
	NTL::ZZ_p::init(NTL::to_ZZ("1234567890123456789012345678901234568123"));
	NTL_ZZ_p::Element a;
	LinBox::integer b, c("123456789012345678901234567890");
	// LinBox::integer b, c("34");
	G.init(a, c );
	G.convert(b, a);
	if ( c != b ) {
		part_pass = false;
                commentator().report (LinBox::Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) << "Error : " << b << " != " << c << std::endl;
	}
        commentator().stop(MSG_STATUS (part_pass));
        pass &= part_pass;

	commentator().stop("NTL_ZZ_p field test suite");
	return pass ? 0 : -1;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

