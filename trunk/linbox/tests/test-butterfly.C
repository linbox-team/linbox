/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-butterfly.C
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
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/util/commentator.h"
#include "linbox/util/vector-factory.h"
#include "linbox/blackbox/butterfly.h"
#include "linbox/switch/boolean.h"
#include "linbox/switch/cekstv.h"

#include "test-common.h"

using namespace LinBox;
using namespace std;

/* Test 1: Boolean switch
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testBooleanSwitch (const Field &F) 
{
	typedef vector <typename Field::Element> Vector;

	commentator.start ("Testing boolean switch", "testBooleanSwitch");

	bool ret = true;

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testBooleanSwitch");

	return ret;
}

/* Test 2: Cekstv switch
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testCekstvSwitch (const Field &F) 
{
	typedef vector <typename Field::Element> Vector;

	commentator.start ("Testing cekstv switch", "testCekstvSwitch");

	bool ret = true;

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testCekstvSwitch");

	return ret;
}

/* Test 3: setButterfly
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testSetButterfly (const Field &F, VectorFactory<Vector> &factory, size_t k) 
{
	typedef std::vector <typename Field::Element> DenseVector;

	commentator.start ("Testing setButterfly", "testSetButterfly", factory.m ());

	bool ret = true, iter_passed;

	Vector v_p;
	DenseVector w (factory.n ()), v1 (factory.n ());
	VectorDomain<Field> VD (F);

	while (factory) {
		commentator.startIteration (factory.j ());

		factory.next (v_p);
		DenseVector v (factory.n ());
		VD.copy (v, v_p);

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector: ";
		VD.write (report, v) << endl;

		std::vector<bool> z (factory.n ());

		for (typename Vector::iterator iter = v_p.begin (); iter != v_p.end (); ++iter)
			z[iter->first] = true;

		BooleanSwitch s (setButterfly (z));
		Butterfly<DenseVector, BooleanSwitch> P (factory.n (), s);

		P.apply (w, v);

		commentator.indent (report);
		report << "Result of apply: ";
		VD.write (report, w) << endl;

		P.applyTranspose (v1, w);

		commentator.indent (report);
		report << "Result of applyTranspose: ";
		VD.write (report, v1) << endl;

		iter_passed = true;

		for (size_t i = 0; i < k; ++i)
			if (F.isZero (w[i]))
				ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Initial block contains zeros" << endl;

		iter_passed = true;

		for (size_t i = k; i < v.size (); ++i)
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

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testSetButterfly");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 1000;
	static integer q = 101U;
	static int iterations = 10;
	static int k = 100;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 1000)", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 101)",   TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",      TYPE_INT,     &iterations },
		{ 'k', "-k K", "K nonzero elements in random vectors (default 100)",   TYPE_INT,     &k },
	};

	typedef Modular<uint16> Field;

	parseArguments (argc, argv, args);
	Field F (q);

	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	cout << "Butterfly preconditioner test suite" << endl << endl;

	RandomSparseMapVectorFactory<Field, NonzeroRandIter<Field> >
		factory (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), n, k, iterations);

	if (!testBooleanSwitch (F)) pass = false;
	if (!testCekstvSwitch  (F)) pass = false;
	if (!testSetButterfly  (F, factory, k)) pass = false;

	return pass ? 0 : -1;
}
