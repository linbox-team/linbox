/* tests/test-modular-float.C
 * Brazenly stolen by bds from
 * tests/test-modular-short.C
 * Brazenly stolen by Zhendong Wan (Copyright (C) 2003) from
 * tests/test-modular.C
 * Copyright (C) 2001, 2002 Bradford Hovinen,
 * Copyright (C) 2002 Dave Saunders
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
 *            Dave Saunders <saunders@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Rename from test-large-modular.C to test-modular.C; made other updates in
 * accordance with changes to Givaro::Modular interace.
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

/*! @file  tests/test-modular-float.C
 * @ingroup tests
 * @brief tests only runFieldTests for modular-float.
 * @test  tests only runFieldTests for modular-float.
 */

#include "givaro/givintprime.h"
#include "linbox/linbox-config.h"
#include "linbox/ring/modular.h"
#include "test-field.h"
using namespace LinBox;

int main (int argc, char **argv)
{
	static size_t n = 10000;
	static unsigned int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	commentator().start("Givaro::Modular<float> field test suite", "Givaro::Modular<float>");
	bool pass = true;

	Givaro::Modular<float> F2 (2);
	Givaro::Modular<float> F5 (5);
	Givaro::Modular<float> F7 (7);
	Givaro::IntPrimeDom IPD;
	integer k = Givaro::Modular<float>::maxCardinality()+1;
	IPD.prevprime(k, k);
	Givaro::Modular<float> Fmax(k);
	IPD.prevprime(k, k/2);
	Givaro::Modular<float> Fmid(k);


	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	std::ostream& report = commentator().report();
	report << "Field F2" << std::endl;
	if (!runFieldTests (F2,  "Givaro::Modular<float>",  iterations, n, false)) pass = false;
	report << "Field F5" << std::endl;
	if (!runFieldTests (F5,  "Givaro::Modular<float>",  iterations, n, false)) pass = false;
	report << "Field F7" << std::endl;
	if (!runFieldTests (F7,  "Givaro::Modular<float>",  iterations, n, false)) pass = false;
	report << "Field Fmax" << std::endl;
	if (!runFieldTests (Fmax,  "Givaro::Modular<float>",  iterations, n, false)) pass = false;
	report << "Field Fmid" << std::endl;
	if (!runFieldTests (Fmid,  "Givaro::Modular<float>",  iterations, n, false)) pass = false;

	commentator().stop("Givaro::Modular<float> field test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
