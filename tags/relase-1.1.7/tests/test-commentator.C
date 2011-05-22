
/* tests/test-commentator.C
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * N.B. This test differs from the others in that the output is not
 * automatically checked for correctness. This is partly because I do not feel
 * it is as critical -- it is easy enough to check oneself, after all -- partly
 * because implementation is rather nontrivial, and partly because of concerns
 * about maintainability. The user must look manually at the output file
 * (specified on the command line) for any anomolies.
 */

#include "linbox/linbox-config.h"

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
	static const char *classes[] =
		{ INTERNAL_DESCRIPTION, INTERNAL_WARNING, INTERNAL_ERROR };
	static const char *levels[] =
		{ "LEVEL_ALWAYS     ", "LEVEL_IMPORTANT  ", "LEVEL_NORMAL     ", "LEVEL_UNIMPORTANT" };

	for (unsigned int i = 0; i < sizeof (classes) / sizeof (char *); ++i)
		for (unsigned int j = 0; j < sizeof (levels) / sizeof (char *); ++j)
			commentator.report (j, classes[i])
				<< "Test string at " << levels[j] << " of class " << classes[i] << endl;
}

// Simple test activity to exercise commentator features

void runTestActivity (bool reportStrings) 
{
	commentator.start ("Test activity", "test", 2);

	if (reportStrings) outputReportStrings ();

	for (int i = 0; i < 2; ++i) {
		commentator.startIteration (i);
		if (reportStrings) outputReportStrings ();

		commentator.start ("Special function 1", "function1");
		if (reportStrings) outputReportStrings ();
		commentator.stop ("done");

		commentator.start ("Special function 2", "function2");
		if (reportStrings) outputReportStrings ();
		commentator.stop ("done");

		commentator.start ("Special function 3", "function3");
		if (reportStrings) outputReportStrings ();
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
	//cout << "Testing primary output...";

	bool ret = true;

	// Disable the brief report for now
	commentator.setMessageClassStream (BRIEF_REPORT, commentator.cnull);

	commentator.setMaxDepth (1);
	runTestActivity (true);

	commentator.setMaxDepth (2);
	commentator.setMaxDetailLevel (Commentator::LEVEL_IMPORTANT);
	runTestActivity (true);

	commentator.setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	runTestActivity (true);

	commentator.setMaxDepth (3);
	runTestActivity (true);

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_ALWAYS);
	runTestActivity (true);

	commentator.setPrintParameters (3, Commentator::LEVEL_IMPORTANT);
	commentator.setPrintParameters (0, Commentator::LEVEL_ALWAYS, "function1");
	commentator.setPrintParameters (10, Commentator::LEVEL_IMPORTANT, "function2");
	commentator.setPrintParameters (10, Commentator::LEVEL_UNIMPORTANT, "function3");
	runTestActivity (true);

	//cout << (ret ? "passed" : "FAILED") << endl;

	return ret;
}

/* Test 2: Brief report
 *
 * Return true on success and false on failure
 */

static bool testBriefReport () 
{
	//cout << "Testing brief report...";

	bool ret = true;

	commentator.setMessageClassStream (BRIEF_REPORT, commentator.report (Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION));
	commentator.setReportStream (commentator.cnull);

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, true, true, true);

	commentator.getMessageClass (BRIEF_REPORT).setMaxDepth (1);
	//cout << "about to ran a testactivity" << endl;
	runTestActivity (false);

	//cout << "ran a testactivity" << endl;

	commentator.getMessageClass (BRIEF_REPORT).setMaxDepth (2);
	runTestActivity (false);

	commentator.getMessageClass (BRIEF_REPORT).setMaxDepth (3);
	runTestActivity (false);

	commentator.setBriefReportParameters (Commentator::OUTPUT_PIPE, true, true, true);

	commentator.getMessageClass (BRIEF_REPORT).setMaxDepth (1);
	runTestActivity (false);

	commentator.getMessageClass (BRIEF_REPORT).setMaxDepth (2);
	runTestActivity (false);

	commentator.getMessageClass (BRIEF_REPORT).setMaxDepth (3);
	runTestActivity (false);

	//cout << (ret ? "passed" : "FAILED") << endl;

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static Argument args[] = {
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.start("Commentator test suite", "commentator");

	if (!testPrimaryOutput ()) pass = false;
	if (!testBriefReport ()) pass = false;

	commentator.stop("commentator test suite");
	//cout << (pass ? "passed" : "FAILED") << endl;
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
