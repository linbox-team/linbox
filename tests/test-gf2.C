/* tests/test-gf2.C
 * Copyright (C) 2003 Bradford Hovinen,
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
 * Updated Mar2016 -bds
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

/*! @file  tests/test-gf2.C
 * @ingroup tests
 * @brief  basic field functionality check
 */

#include "linbox/linbox-config.h"
#include "linbox/field/gf2.h"
#include "test-field.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static unsigned int n = 100;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to NxN.",      TYPE_INT,     &n },
		END_OF_ARGUMENTS
	};
	parseArguments (argc, argv, args);

	bool pass = true;
	commentator().start("GF2 field test suite", "GF2");

	GF2 F;

	pass = pass and testField (F, "GF2");
	pass = pass and runFieldTests (F, "GF2", 1, n, false);

	commentator().stop("GF2 field test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
