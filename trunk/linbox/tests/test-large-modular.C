/* -*- mode: c; style: linux -*- */

/* tests/test-diagonal.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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
	static size_t n = 10;
	static integer q = 4294967291U;
	static int iterations = 100;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)",        TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 4294967291)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",          TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);

	srand (time (NULL));

	cout << "LargeModular field test suite" << endl << endl;
	bool res, pass = true;

	LargeModular F (q);

	test_header("LargeModular field", report);
	res = test_field<LargeModular>    (F, n, report, iterations);
	test_trailer(res, report);
	pass &= res;

	Field_archetype K(new LargeModular(101));
	test_header("Arch of Env of LargeModular field", report);
	res = test_field<Field_archetype>    (K, n, report, iterations);
	test_trailer(pass, report);
	pass &= res;

	return pass ? 0 : -1;
}
