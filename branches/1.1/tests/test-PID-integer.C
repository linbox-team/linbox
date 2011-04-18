/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* tests/test-PID-integer.C
 * Copyright (C) 2011 Dave Saunders
 *
 * morphed from another test by Dave Saunders <saunders@cis.udel.edu>
 *
 * See COPYING for license information.
 */

/*! @file  tests/test-PID-integer.C
 * @ingroup tests
 * @brief run runFieldTests testRandomIterator tests on PID_integer
 * @test run runFieldTests testRandomIterator tests on PID_integer
 */


#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>

#include "test-common.h"
#include "test-field.h"

#include "linbox/field/PID-integer.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static size_t n = 10000;
	static int iterations = 1;
	static int trials = 100000;
	static int categories = 100;
	static int hist_level = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	commentator.start("PID-integer ring test suite", "PID-integer");
	bool pass = true;

	PID_integer ZZ;

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	if (!runBasicRingTests (ZZ,  "PID-integer",  iterations, false)) pass = false;
//	if (!testRandomIterator (ZZ,  "PID-integer", trials, categories, hist_level)) pass = false;

	commentator.stop("PID-integer, field test suite");
	return pass ? 0 : -1;
}
