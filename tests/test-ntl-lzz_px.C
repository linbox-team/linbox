/* tests/test-ntl-lzz_pX.C
 * evolved by bds from tests/test-ntl-lzz_pE.C
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

/*! @file  tests/test-ntl-lzz_pX.C
 * @ingroup tests
 * @brief no doc.
 * @test no doc.
 */


#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "linbox/ring/ntl.h"

#include "test-field.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static integer q = 31;
	static size_t n = 10000;
	static int iterations = 1;

	static Argument args[] = {
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	commentator().start("NTL_zz_pX field test suite", "NTL_zz_p");
	bool pass = true;

	NTL_zz_pX R(q);

    integer ch;
	if (R.characteristic(ch) != q)
		exit(-1);

	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	if (!runBasicRingTests (R, "NTL_zz_pX", (unsigned int)iterations, false)) pass = false;
	if (!runPIRTests (R, "NTL_zz_pX", (unsigned int)iterations, false)) pass = false;

	// Test generic i/o

    std::stringstream ss;
    ss << 3 << ' ' << Givaro::Integer(1) << ' ' << Givaro::Integer(-2)
       << ' ' << Givaro::Integer(5) << ' ' << Givaro::Integer(-7) << std::endl;

    NTL_zz_pX::Element P;

    R.read(ss, P);
    R.write(std::clog << "# P: ", P) << std::endl;


	// needs PID tests as well...

#if 0
	FieldArchetype K(new NTL_zz_pX(101));

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of UnField<NTL::zz_pX> field"))
		pass = false;
#endif

	commentator().stop("NTL_zz_pX field test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
