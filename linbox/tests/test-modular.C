/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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
 * See COPYING for license information.
 */

#include "linbox-config.h"

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
	static size_t n = 10000;
	static int iterations = 10;
	static int trials = 1000000;
	static int categories = 100;
	static int hist_level = 1;

	static Argument args[] = {
		{ 'K', "-K Q", "Operate over the \"field\" GF(Q) [1] for integer modulus (default 18446744073709551557)", TYPE_INTEGER, &q1 },
		{ 'Q', "-Q Q", "Operate over the \"field\" GF(Q) [1] for uint32 modulus (default 2147483647)", TYPE_INTEGER, &q2 },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint16 modulus (default 65521)", TYPE_INTEGER, &q3 },
		{ 'p', "-p P", "Operate over the \"field\" GF(Q) [1] for uint8 modulus (default 101)", TYPE_INT, &q4 },
		{ 'n', "-n N", "Set dimension of test vectors to NxN (default 10000)",      TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",           TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test (default 1000000)", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test (default 100)", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test (default 1)", TYPE_INT, &hist_level },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	cout << endl << "Modular field test suite" << endl;
	cout.flush ();
	bool pass = true;

	Modular<integer> F_integer (q1);
	Modular<uint32> F_uint32 ((uint32) q2);
	Modular<uint16> F_uint16 ((uint16) q3);
	Modular<uint8> F_uint8 ((uint8) q4);
	Modular<float> F_float ((float) q4);

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	if (!runFieldTests (F_integer, "Modular<integer>", iterations, n, false)) pass = false;
	if (!runFieldTests (F_uint32,  "Modular<uint32>",  iterations, n, false)) pass = false;
	if (!runFieldTests (F_uint16,  "Modular<uint16>",  iterations, n, false)) pass = false;
	if (!runFieldTests (F_uint8,  "Modular<uint8>",  iterations, n, false)) pass = false;
	if (!runFieldTests (F_float,  "Modular<float>",  iterations, n, false)) pass = false;

	//if (!testRandomIterator (F_integer, "Modular<integer>", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F_uint32,  "Modular<uint32>", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F_uint16,  "Modular<uint16>", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F_uint8,  "Modular<uint8>", trials, categories, hist_level)) pass = false;

	return pass ? 0 : -1;
}
