/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-commentator.C
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 *
 * See COPYING for license information.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include "linbox/util/commentator.h"

#include "test-common.h"

using namespace LinBox;
using namespace std;

// Output various report strings

void outputReportStrings () 
{
	static char *classes[] =
		{ INTERNAL_DESCRIPTION, INTERNAL_WARNING, INTERNAL_ERROR };
	static char *levels[] =
		{ "LEVEL_ALWAYS     ", "LEVEL_IMPORTANT  ", "LEVEL_NORMAL     ", "LEVEL_UNIMPORTANT" };

	for (unsigned int i = 0; i < sizeof (classes) / sizeof (char *); ++i)
		for (unsigned int j = 0; j < sizeof (levels) / sizeof (char *); ++j)
			commentator.report (j + 1, classes[i])
				<< "Test string at " << levels[j] << " of class " << classes[i] << endl;
}

// Simple test activity to exercise commentator features

void runTestActivity () 
{
	commentator.start ("Test activity", "test", 2);

	outputReportStrings ();

	for (int i = 0; i < 2; ++i) {
		commentator.startIteration (i);
		outputReportStrings ();

		commentator.start ("Special function", "function1");
		outputReportStrings ();
		commentator.stop ("done");

		commentator.start ("Another special function", "function2");
		outputReportStrings ();
		commentator.stop ("done");

		commentator.start ("Another special function", "function3");
		outputReportStrings ();
		commentator.stop ("done");

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop ("done", (const char *) 0, "test");
}

/* Test 1: Primary output
 *
 * Return true on success and false on failure
 */

static bool testPrimaryOutput () 
{
	cout << "Testing primary output...";

	bool ret = true;

	// Disable the brief report for now
	commentator.setMessageClassStream (BRIEF_REPORT, commentator.cnull);

	commentator.setMaxDepth (1);
	runTestActivity ();

	commentator.setMaxDepth (2);
	commentator.setMaxDetailLevel (Commentator::LEVEL_IMPORTANT);
	runTestActivity ();

	commentator.setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	runTestActivity ();

	commentator.setMaxDepth (3);
	runTestActivity ();

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_ALWAYS);
	runTestActivity ();

	commentator.setPrintParameters (3, Commentator::LEVEL_IMPORTANT);
	commentator.setPrintParameters (0, Commentator::LEVEL_ALWAYS, "function1");
	commentator.setPrintParameters (10, Commentator::LEVEL_IMPORTANT, "function2");
	commentator.setPrintParameters (10, Commentator::LEVEL_UNIMPORTANT, "function3");
	runTestActivity ();

	cout << (ret ? "passed" : "FAILED") << endl;

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static Argument args[] = {
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	cout << "Commentator test suite" << endl << endl;

	if (!testPrimaryOutput ()) pass = false;

	return pass ? 0 : -1;
}
