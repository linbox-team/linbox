/* -*- mode: c; style: linux -*- */

/* tests/test-diagonal.C
 * Copyright (C) 2001, 2002 Bradford Hovinen,
 * Copyright (C) 2002 Dave Saunders
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
 *            Dave Saunders <saunders@cis.udel.edu>
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/large-modular.h"

#include "test-common.h"
#include "test-field-common.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	ofstream report;
	static integer q = 4294967291U;

	static Argument args[] = {
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 4294967291)", TYPE_INTEGER, &q },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	cout << "LargeModular field test suite" << endl << endl;
	cout.flush ();
	bool pass = true;

	LargeModular F (q);

	if (!testField<LargeModular> (F, "Testing LargeModular field"))
		pass = false;

#if 0
	Field_archetype K(new LargeModular(101));

	if (!testField<Field_archetype> (K, "Testing archetype with envelope of LargeModular field"))
		pass = false;
#endif

	return pass ? 0 : -1;
}
