/* -*- mode: c; style: linux -*- */

/* tests/test-unparametric-field.cpp
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/unparametric.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static integer q = 4294967291U;

	static Argument args[] = {
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	cout << "Unparametrix<double> field test suite" << endl << endl;
	cout.flush ();
	bool pass = true;

	UnparametricField<double> F;

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	if (!testField<UnparametricField<double> > (F, "Testing UnparametricField<double> field"))
		pass = false;

#if 0
	FieldArchetype K(new UnparametricField(101));

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of UnField field"))
		pass = false;
#endif

	return pass ? 0 : -1;
}
