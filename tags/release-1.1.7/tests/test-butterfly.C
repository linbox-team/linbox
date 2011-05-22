
/* tests/test-butterfly.C
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 *
 * See COPYING for license information.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/randiter/nonzero.h"
#include "linbox/util/commentator.h"
#include "linbox/vector/stream.h"
#include "linbox/blackbox/butterfly.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/submatrix.h"
#include "linbox/solutions/det.h"
#include "linbox/switch/boolean.h"
#include "linbox/switch/cekstv.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;
using namespace std;

/* Test 1: setButterfly/BooleanSwitch
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testSetButterfly (const Field &F, VectorStream<Vector> &stream, size_t k) 
{
	commentator.start ("Testing setButterfly", "testSetButterfly", stream.size ());

	bool ret = true, iter_passed;

	unsigned long real_k;

	Vector v_p;
	typename LinBox::Vector<Field>::Dense w (stream.dim ()), v1 (stream.dim ());
	VectorDomain<Field> VD (F);

	while (stream) {
		commentator.startIteration (stream.j ());

		stream >> v_p;
		typename LinBox::Vector<Field>::Dense v (stream.n ());
		VD.copy (v, v_p);

		real_k = v_p.first.size ();

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector: ";
		VD.write (report, v) << endl;

		std::vector<bool> z (stream.n ());

		for (typename Vector::first_type::iterator iter = v_p.first.begin (); iter != v_p.first.end (); ++iter)
			z[*iter] = true;

		vector<bool> z_vec = setButterfly (z);

		BooleanSwitchFactory factory (z_vec);
		Butterfly<Field, BooleanSwitch> P (F, stream.dim (), factory);

		P.apply (w, v);

		report << "Result of apply: ";
		VD.write (report, w) << endl;

		P.applyTranspose (v1, w);

		report << "Result of applyTranspose: ";
		VD.write (report, v1) << endl;

		iter_passed = true;

		for (size_t i = 0; i < real_k; ++i)
			if (F.isZero (w[i]))
				ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Initial block contains zeros" << endl;

		iter_passed = true;

		for (size_t i = real_k; i < v.size (); ++i)
			if (!F.isZero (w[i]))
				ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Nonzero entries outside initial block" << endl;

		if (!VD.areEqual (v, v1)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: P^T != P^-1" << endl;
			ret = false;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	stream.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testSetButterfly");

	return ret;
}

/* Test 2: Cekstv switch
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testCekstvSwitch (const Field &F, unsigned int iterations, size_t n, size_t r) 
{
	commentator.start ("Testing cekstv switch", "testCekstvSwitch", iterations);

	bool ret = true;

	unsigned int failures = 0;

	typename Vector<Field>::SparsePar d1;

	RandomSparseStream<Field, typename Vector<Field>::SparsePar> stream (F, (double) r / (double) n, n, iterations);
	VectorDomain<Field> VD (F);

	unsigned long real_r;
	typename Field::Element det_Ap;

	typename Field::Element one;

	F.init (one, 1);

	while (stream) {
		commentator.startIteration (stream.pos ());

		stream >> d1;

		typename Vector<Field>::Dense d (n);
		VD.copy (d, d1);

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector: ";
		VD.write (report, d) << endl;

		real_r = d1.first.size ();

		report << "Real rank: " << real_r << endl;

		typename Field::RandIter r (F);
		typename CekstvSwitch<Field>::Factory factory (r);
		Butterfly<Field, CekstvSwitch<Field> > P (F, n, factory);
		Butterfly<Field, CekstvSwitch<Field> > Q (F, n, factory);
		typedef Butterfly<Field, CekstvSwitch<Field> > Blackbox1;

		Diagonal<Field> D (F, d);
		typedef Diagonal<Field> Blackbox2;

		Compose<Blackbox1>  DQ (&P, &Q);
		typedef Compose<Blackbox1, Blackbox1> Blackbox3;
		Compose<Blackbox1, Blackbox3> A (P, DQ);
		typedef Compose<Blackbox1, Blackbox3> Blackbox4;
		//Compose<typename Vector<Field>::Dense> A (&P, &DQ);

		Submatrix<Blackbox4> Ap (&A, 0, 0, real_r, real_r);

		det (det_Ap, Ap,  Method::Wiedemann());

		report << "Deteriminant of r x r leading principal minor: ";
		F.write (report, det_Ap) << endl;

		if (F.isZero (det_Ap)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
				<< "WARNING: Determinant of r x r leading principal minor is zero" << endl;
			++failures;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
		<< "Total failures: " << failures << endl;

	// FIXME: I need a theoretical error bound
	if (failures > iterations / 5) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Too many failures. This is likely a bug." << endl;
		ret = false;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testCekstvSwitch");

	return ret;
}

/* Test 3: Random linearity
 *
 * Construct a random dense matrix and a submatrix thereof. Call testLinearity
 * in test-generic.h to test that the submatrix is a linear operator
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrices
 * iterations - Number of iterations to run
 * N - Number of random vectors to which to apply
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testRandomLinearity (const Field                                 &F,
				 VectorStream<typename Vector<Field>::Dense> &v1_stream,
				 VectorStream<typename Vector<Field>::Dense> &v2_stream) 
{
	commentator.start ("Testing random linearity", "testRandomLinearity", v1_stream.size ());

	typename Field::RandIter r (F);
	typename CekstvSwitch<Field>::Factory factory (r);
	Butterfly<Field, CekstvSwitch<Field> > A (F, v1_stream.dim (), factory);

	bool ret = testLinearity (F, A, v1_stream, v2_stream);

	v1_stream.reset ();
	v2_stream.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomLinearity");

	return ret;
}

/* Test 4: Random transpose
 *
 * Construct a random dense matrix and a submatrix thereof. Call testLinearity
 * in test-generic.h to test that the submatrix is a linear operator
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrices
 * iterations - Number of iterations to run
 * N - Number of random vectors to which to apply
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testRandomTranspose (const Field                                 &F,
				 VectorStream<typename Vector<Field>::Dense> &v1_stream,
				 VectorStream<typename Vector<Field>::Dense> &v2_stream) 
{
	commentator.start ("Testing random transpose", "testRandomTranspose", v1_stream.size ());

	typename Field::RandIter r (F);
	typename CekstvSwitch<Field>::Factory factory (r);
	Butterfly<Field, CekstvSwitch<Field> > A (F, v1_stream.dim (), factory);

	bool ret = testTranspose (F, A, v1_stream, v2_stream);

	v1_stream.reset ();
	v2_stream.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomTranspose");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 100;
	static integer q = 2147483647U;
	static int iterations = 1; // was 10
	static int k = 10;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.",      TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.",           TYPE_INT,     &iterations },
		{ 'k', "-k K", "K nonzero elements in random vectors.",        TYPE_INT,     &k },
    	{ '\0' }
	};

	typedef Modular<LinBox::uint32> Field;

	parseArguments (argc, argv, args);
	Field F (q);

	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	commentator.start("Butterfly preconditioner test suite", "butterfly preconditioner");

	RandomSparseStream<Field, Vector<Field>::Sparse, NonzeroRandIter<Field> >
		stream (F, NonzeroRandIter<Field> (F, Field::RandIter (F)),
			(double) k / (double) n, n, iterations);
	RandomDenseStream<Field> v1_stream (F, n, iterations);
	RandomDenseStream<Field> v2_stream (F, n, iterations);

	if (!testSetButterfly  (F, stream, k)) pass = false;
	if (!testCekstvSwitch  (F, iterations, n, k)) pass = false;
	if (!testRandomLinearity (F, v1_stream, v2_stream)) pass = false;
	if (!testRandomTranspose (F, v1_stream, v2_stream)) pass = false;

	commentator.stop("butterfly preconditioner test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
