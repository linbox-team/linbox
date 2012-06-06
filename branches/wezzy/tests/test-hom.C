/* tests/test-hom.C
 * Copyright (C) LinBox
 *
 * Written by Dave Saunders <saunders@cis.udel.edu>
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


/*! @file  tests/test-hom.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */



#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include "linbox/field/modular.h"
#include "linbox/field/hom.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static integer q = 65521U;
	static size_t n = 10000;
	static int iterations = 10;
	static int trials = 100000;
	static int categories = 100;
	static int hist_level = 1;

	static Argument args[] = {
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint16_t modulus.", TYPE_INTEGER, &q },
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	commentator().start("Hom test suite", "Hom");
	bool pass = true;

	Modular<uint32_t> F_uint32_t ((uint32_t) q);
	Modular<uint16_t> F_uint16_t ((uint16_t) q);
	Hom<Modular<uint16_t>, Modular<uint32_t> > iso(F_uint16_t, F_uint32_t);

	uint16_t x=2, y;
	uint32_t z=2, w;
	iso.image(w, x);
	pass = pass && F_uint32_t.areEqual(z, w);
	iso.preimage(y, z);
	pass = pass && F_uint16_t.areEqual(x, y);

	/* for image field!
	uint32_t x, y, z, w;
	iso.smul(x, 2, 3);
	iso.mul(y, 2, 3);
	pass = pass && iso.areEqual(x, y);
	iso.smulin(x, 5);
	iso.mulin(y, 5);
	pass = pass && iso.areEqual(x, y);

	iso.saxpy(z, 7, x, 11);
	iso.axpy(w, 7, x, 11);
	pass = pass && iso.areEqual(z, w);
	iso.saxpyin(z, 7, x);
	iso.axpyin(w, 7, x);
	pass = pass && iso.areEqual(z, w);
	*/

	commentator().stop("Hom test suite");
	return pass ? 0 : -1;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

