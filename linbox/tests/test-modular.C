/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* tests/test-modular.C
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
 * accordance with changes to Modular interace.
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
 * License along with this library; if not, write to the Free Software Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */


/*! @file  tests/test-modular.C
 * @ingroup tests
 *
 * @brief  runFieldTests on various <code>Modular<XXX></code> fields.
 * @test   runFieldTests on various <code>Modular<XXX></code> fields.
 */



#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>

#include "linbox/field/modular.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{
	static integer q1("18446744073709551557");
	static integer q2 = 2147483647U;
	static integer q3 = 65521U;
	static int q4 = 101;
	static size_t n = 500;
	static int iterations = 1;
	static int trials = 100000;
	static int categories = 100;
	static int hist_level = 1;

	static Argument args[] = {
		{ 'K', "-K Q", "Operate over the \"field\" GF(Q) [1] for integer modulus.", TYPE_INTEGER, &q1 },
		{ 'Q', "-Q Q", "Operate over the \"field\" GF(Q) [1] for uint32_t modulus.", TYPE_INTEGER, &q2 },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint16_t modulus.", TYPE_INTEGER, &q3 },
		{ 'p', "-p P", "Operate over the \"field\" GF(Q) [1] for uint8_t modulus.", TYPE_INT, &q4 },
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	commentator.start("Modular test suite", "Modular ");
	bool pass = true;

	Modular<integer> F_integer (q1);
	Modular<uint32_t> F_uint32_t ((uint32_t) q2);
	Modular<uint16_t> F_uint16_t ((uint16_t) q3);
	Modular<uint8_t> F_uint8_t ((uint8_t) q4);
	Modular<float> F_float ((float) q4);


	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	if (!runFieldTests (F_integer,"Modular<integer>", iterations, n, false)) pass = false;
	if (!runFieldTests (F_uint32_t, "Modular<uint32_t>",  iterations, n, false)) pass = false;
	if (!runFieldTests (F_uint16_t, "Modular<uint16_t>",  iterations, n, false)) pass = false;
	if (!runFieldTests (F_uint8_t,  "Modular<uint8_t>",   iterations, n, false)) pass = false;
	if (!runFieldTests (F_float,  "Modular<float>",   iterations, n, false)) pass = false;

	if (!testRandomIterator (F_integer, "Modular<integer>", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F_uint32_t,  "Modular<uint32_t>",  trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F_uint16_t,  "Modular<uint16_t>",  trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F_uint8_t,    "Modular<uint8_t>",  trials, categories, hist_level)) pass = false;

	commentator.stop(MSG_STATUS (pass), "Modular test suite");
	return pass ? 0 : -1;
}
