/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-gf2.C
 * Copyright (C) 2003 Bradford Hovinen,
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
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

#include "linbox/field/gf2.h"
#include "linbox/field/modular.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

#if 0

template <class Vector1, class Vector2>
static bool testDotProduct (const GF2 &F, LinBox::VectorStream<Vector1> &stream1,
			    LinBox::VectorStream<Vector2> &stream2) 
{
	LinBox::commentator.start ("Testing GF2 dot product", "testDotProduct", stream1.size ());

	bool ret = true;

	Vector1 v1;
	Vector2 v2;

	Modular<uint16> MF (2);
	VectorDomain<Modular<uint16> > MF_VD (MF);

	LinBox::Vector<Modular<uint16> >::Dense v3 (stream1.dim ()), v4 (stream1.dim ());

	typename Field::Element sigma, rho;

	LinBox::VectorDomain<Field> VD (F);

	size_t j;

	LinBox::VectorWrapper::ensureDim (v1, stream1.n ());
	LinBox::VectorWrapper::ensureDim (v2, stream2.n ());

	LinBox::Timer timer;
	double totaltime = 0.0;

	while (stream1 && stream2) {
		LinBox::commentator.startIteration (stream1.j ());

		F.init (sigma, 0);

		stream1.next (v1);
		stream2.next (v2);

		for (j = 0; j < stream1.n (); j++)
			F.axpyin (sigma,
				  LinBox::VectorWrapper::constRef<Field> (v1, j),
				  LinBox::VectorWrapper::constRef<Field> (v2, j));

		std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1:  ";
		VD.write (report, v1) << std::endl;

		report << "Input vector 2:  ";
		VD.write (report, v2) << std::endl;

		timer.start ();
		VD.dot (rho, v1, v2);
		timer.stop ();
		totaltime += timer.realtime ();

		report << "True dot product: ";
		F.write (report, sigma) << std::endl;

		report << "Dot product from vector domain: ";
		F.write (report, rho) << std::endl;

		if (!F.areEqual (sigma, rho)) {
			ret = false;
			LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Dot products are not equal" << std::endl;
		}

		LinBox::commentator.stop ("done");
		LinBox::commentator.progress ();
	}

	LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
		<< "Average time for dot product: " << totaltime / stream1.m () << std::endl;

	LinBox::commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDotProduct");

	stream1.reset ();
	stream2.reset ();

	return ret;
}

#endif

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

	cout << "GF2 field test suite" << endl << endl;
	cout.flush ();
	bool pass = true;

	GF2 F;

	RandomDenseStream<GF2> stream1 (F, n, iterations), stream2 (F, n, iterations);
	RandomSparseStream<GF2> stream3 (F, 0.1, n, iterations), stream4 (F, 0.1, n, iterations);

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	if (!runFieldTests (F, "GF2", iterations, n, false)) pass = false;
	if (!testRandomIterator (F, "GF2", trials, categories, hist_level)) pass = false;

	commentator.start ("Testing VectorDomain <GF2>", "main");

#if 0
	if (!testDotProduct (F, "dense/dense", stream1, stream2)) pass = false;
	if (!testDotProduct (F, "dense/sparse", stream1, stream3)) pass = false;
	if (!testDotProduct (F, "sparse/dense", stream3, stream1)) pass = false;
	if (!testDotProduct (F, "sparse/sparse", stream3, stream4)) pass = false;
#endif

	if (!testAddMul (F, "dense", stream1, stream2)) pass = false;
	if (!testAddMul (F, "sparse", stream3, stream4)) pass = false;

	if (!testSubMul (F, "dense", stream1, stream2)) pass = false;
	if (!testSubMul (F, "sparse", stream3, stream4)) pass = false;

	if (!testAXPY (F, "dense", stream1, stream2)) pass = false;
	if (!testAXPY (F, "sparse", stream3, stream4)) pass = false;

	if (!testCopyEqual (F, "dense/dense", stream1, stream2)) pass = false;
	if (!testCopyEqual (F, "dense/sparse", stream1, stream3)) pass = false;
	if (!testCopyEqual (F, "sparse/dense", stream3, stream1)) pass = false;
	if (!testCopyEqual (F, "sparse/sparse", stream3, stream4)) pass = false;

	commentator.stop (MSG_STATUS (pass), (const char *) 0, "main");

#if 0
	FieldArchetype K(new LargeModular(101));

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of Modular field"))
		pass = false;
#endif

	return pass ? 0 : -1;
}
