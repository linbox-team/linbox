/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-modular-double.C
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

#include "linbox/field/modular-double.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static size_t n = 10000;
	static int iterations = 10;
	static int trials = 1000000;
	static int categories = 100;
	static int hist_level = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to NxN (default 10000)",      TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",           TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test (default 1000000)", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test (default 100)", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test (default 1)", TYPE_INT, &hist_level },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	cout << endl << "Modular<double> field test suite" << endl;
	cout.flush ();
	bool pass = true;

	//Modular<double> F_int (32749); // largest
	//Modular<double> F2_int (2); // largest
	Modular<double> F3_int (3); // largest
	Modular<double> F5_int (5); // largest
	Modular<double> F7_int (7); // largest
	//Modular<double> F11_int (11); // largest
	//Modular<double> G_int (65521); // largest
	//Modular<double> H_int (1099511627689); // largest

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	//if (!runFieldTests (F2_int,  "Modular<double>",  iterations, n, false)) pass = false;
	if (!runFieldTests (F3_int,  "Modular<double>",  iterations, n, false)) pass = false;
	if (!runFieldTests (F5_int,  "Modular<double>",  iterations, n, false)) pass = false;
	if (!runFieldTests (F7_int,  "Modular<double>",  iterations, n, false)) pass = false;
	//if (!runFieldTests (F11_int,  "Modular<double>",  iterations, n, false)) pass = false;
	//if (!testRandomIterator (F_int,  "Modular<double>", trials, categories, hist_level)) pass = false;
	//if (!runFieldTests (G_int,  "Modular<double>",  iterations, n, false)) pass = false;
	//if (!testRandomIterator (G_int,  "Modular<double>", trials, categories, hist_level)) pass = false;
	//if (!runFieldTests (H_int,  "Modular<double>",  iterations, n, false)) pass = false;
	//if (!testRandomIterator (H_int,  "Modular<double>", trials, categories, hist_level)) pass = false;

	return pass ? 0 : -1;
}
