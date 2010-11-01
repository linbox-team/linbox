
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

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>

#include "linbox/field/gf2.h"
#include "linbox/field/modular.h"

#include "test-generic-for-quad.h"
// #include "test-generic.h"
// #include "test-vector-domain.h"

using namespace LinBox;
using LinBox::uint16;
using LinBox::uint32;

static bool testDotProductGF2 (const GF2 &F, const char *, //desc,
			       VectorStream<Vector<GF2>::Dense> &stream1,
			       VectorStream<Vector<GF2>::Dense> &stream2) 
{
	LinBox::commentator.start ("Testing GF2 dot product (dense/dense)", "testDotProduct", stream1.size ());

	bool ret = true;

	Vector<GF2>::Dense v1, v2;

	Modular<uint16> MF (2);
	VectorDomain<Modular<uint16> > MF_VD (MF);

	LinBox::Vector<Modular<uint16> >::Dense w1 (stream1.dim ()), w2 (stream1.dim ());

	uint16 sigma;
	bool rho;

	LinBox::VectorDomain<GF2> VD (F);

	v1.resize (stream1.dim ());
	v2.resize (stream2.dim ());

	Vector<GF2>::Dense::const_iterator i1;
	Vector<GF2>::Dense::const_iterator i2;
	Vector<Modular<uint16> >::Dense::iterator j1, j2;

	Timer timer;
	double totaltime = 0.0;

	while (stream1 && stream2) {
		LinBox::commentator.startIteration (stream1.j ());

		stream1.next (v1);
		stream2.next (v2);

		// Copy v1 and v2 into w1 and w2
		for (i1 = v1.begin (), j1 = w1.begin (); i1 != v1.end (); ++i1, ++j1) {
			if (*i1)
				*j1 = 1;
			else
				*j1 = 0;
		}

		for (i2 = v2.begin (), j2 = w2.begin (); i2 != v2.end (); ++i2, ++j2) {
			if (*i2)
				*j2 = 1;
			else
				*j2 = 0;
		}

		std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1:  ";
		VD.write (report, v1) << std::endl;

		report << "Input vector 2:  ";
		VD.write (report, v2) << std::endl;

		timer.start ();
		VD.dot (rho, v1, v2);
		timer.stop ();
		totaltime += timer.realtime ();

		MF_VD.dot (sigma, w1, w2);
		sigma &= 1;

		report << "True dot product: ";
		F.write (report, sigma) << std::endl;

		report << "Dot product from vector domain: ";
		F.write (report, rho) << std::endl;

		if ((sigma && !rho) || (rho && !sigma)) {
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

static bool testDotProductGF2 (const GF2 &F, const char *, //desc,
			       VectorStream<Vector<GF2>::Dense> &stream1,
			       VectorStream<Vector<GF2>::Sparse> &stream2) 
{
	LinBox::commentator.start ("Testing GF2 dot product (dense/sparse)", "testDotProduct", stream1.size ());

	bool ret = true;

	Vector<GF2>::Dense v1;
	Vector<GF2>::Sparse v2;

	Modular<uint16> MF (2);
	VectorDomain<Modular<uint16> > MF_VD (MF);

	LinBox::Vector<Modular<uint16> >::Dense w1 (stream1.dim ());
	LinBox::Vector<Modular<uint16> >::SparseSeq w2;

	uint16 sigma;
	bool rho;

	LinBox::VectorDomain<GF2> VD (F);

	v1.resize (stream1.dim ());

	Vector<GF2>::Dense::const_iterator i1;
	Vector<GF2>::Sparse::const_iterator i2;
	Vector<Modular<uint16> >::Dense::iterator j1;

	Timer timer;
	double totaltime = 0.0;

	while (stream1 && stream2) {
		LinBox::commentator.startIteration (stream1.j ());

		stream1.next (v1);
		stream2.next (v2);

		// Copy v1 and v2 into w1 and w2
		for (i1 = v1.begin (), j1 = w1.begin (); i1 != v1.end (); ++i1, ++j1) {
			if (*i1)
				*j1 = 1;
			else
				*j1 = 0;
		}

		w2.clear ();

		for (i2 = v2.begin (); i2 != v2.end (); ++i2)
			w2.push_back (std::pair<size_t, uint16> (*i2, 1));

		std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1:  ";
		VD.write (report, v1) << std::endl;

		report << "Input vector 2:  ";
		VD.write (report, v2) << std::endl;

		timer.start ();
		VD.dot (rho, v1, v2);
		timer.stop ();
		totaltime += timer.realtime ();

		MF_VD.dot (sigma, w1, w2);
		sigma &= 1;

		report << "True dot product: ";
		F.write (report, sigma) << std::endl;

		report << "Dot product from vector domain: ";
		F.write (report, rho) << std::endl;

		if ((sigma && !rho) || (rho && !sigma)) {
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

int main (int argc, char **argv)
{
	static unsigned int n = 1000;
	static int iterations = 2;
	static int trials = 100000;
	static int categories = 100;
	static int hist_level = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to NxN.",      TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.",           TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.start("GF2 field test suite", "GF2");
	bool pass = true;

	GF2 F;

	uint32 seed = time (NULL);

	RandomDenseStreamGF2 stream1 (F, seed, n, iterations), stream2 (F, seed ^ 0xdeadbeef, n, iterations);
	RandomSparseStreamGF2<Vector<GF2>::Sparse>
		stream3 (F, seed + 2, 0.1, n, iterations),
		stream4 (F, seed + 3, 0.1, n, iterations);

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	commentator.start ("Testing GF2", "main", 10);
	
	if (!testFieldNegation         (F, "GF2", iterations))      pass = false; commentator.progress ();
	if (!testFieldInversion        (F, "GF2", iterations))      pass = false; commentator.progress ();
	if (!testFieldAxioms           (F, "GF2", iterations))      pass = false; commentator.progress ();
	if (!testFieldAssociativity    (F, "GF2", iterations))      pass = false; commentator.progress ();
	if (!testFieldCharacteristic   (F, "GF2", iterations))      pass = false; commentator.progress ();
	if (!testGeometricSummation    (F, "GF2", iterations, 100)) pass = false; commentator.progress ();
	if (!testFreshmansDream        (F, "GF2", iterations))      pass = false; commentator.progress ();
	if (!testArithmeticConsistency (F, "GF2", iterations))      pass = false; commentator.progress ();
	if (!testAxpyConsistency       (F, "GF2", iterations))      pass = false; commentator.progress ();

	commentator.stop (MSG_STATUS (pass), (const char *) 0, "main");

	if (!testRandomIterator (F, "GF2", trials, categories, hist_level)) pass = false;

	commentator.start ("Testing VectorDomain <GF2>", "main");

	if (!testDotProductGF2 (F, "dense/dense", stream1, stream2)) pass = false;
	if (!testDotProductGF2 (F, "dense/sparse", stream1, stream3)) pass = false;

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

	commentator.stop("GF2 field test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
