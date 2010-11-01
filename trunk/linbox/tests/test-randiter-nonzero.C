
/* tests/test-randiter-nonzero.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/randiter/modular.h"
#include "linbox/randiter/nonzero.h"

#include "test-common.h"

using namespace std;
using namespace LinBox;

/* Test 1: Nonzero random elements
 *
 * Construct a sequence of random elements from the nonzero random iterator and
 * check that they are all nonzero
 *
 * F - Field over which to perform computations
 * iterations - Number of iterations to perform
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testNonzeroRandom (Field &F, int iterations) 
{
	int i;

	commentator.start ("Testing nonzero random elements", "testNonzeroRandom", iterations);

	bool ret = true;

	typename Field::Element x;
	typename Field::RandIter r (F);

	NonzeroRandIter <Field, typename Field::RandIter> rp (F, r);

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		rp.random (x);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Element: ";
		F.write (report, x);
		report << endl;

		if (F.isZero (x)) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Element is zero" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testNonzeroRandom");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static integer q = 101;
	static int iterations = 1000;

	static Argument args[] = {
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);
	Modular<uint32> F (q);

	srand (time (NULL));

	commentator.start("Nonzero random iterator test suite", "NonzeroRandIter");

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	if (!testNonzeroRandom (F, iterations)) pass = false;

	commentator.stop("Nonzero random iterator test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
