/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-sparse.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/vector-factory.h"
#include "linbox/field/vector-domain.h"
#include "linbox/randiter/nonzero.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

/* Test 1: Application of identity matrix onto random vectors
 *
 * Construct the identity matrix and a series of randomly-generated
 * vectors. Apply the identity to each vector and test whether the input and
 * output are equal.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 *
 * Return true on success and false on failure
 */

template <class Row, class Field, class Vector>
static bool testIdentityApply (Field &F, const char *text, VectorFactory<Vector> &factory) 
{
	typedef SparseMatrix0 <Field, Vector, Row> Blackbox;

	ostringstream str;
	str << "Testing identity apply (" << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testIdentityApply", factory.m ());

	bool ret = true;
	bool iter_passed = true;

	VectorDomain<Field> VD (F);
	StandardBasisFactory<Field, Row> f1 (F, factory.n ());
	Blackbox A (F, f1);

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Matrix:" << endl;
	A.write (report, Blackbox::FORMAT_PRETTY);

	Vector v, w;

	VectorWrapper::ensureDim (v, factory.n ());
	VectorWrapper::ensureDim (w, factory.n ());

	while (factory) {
		commentator.startIteration (factory.j ());

		iter_passed = true;

		factory.next (v);

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		VD.write (report, v);
		report << endl;

		A.apply (w, v);

		commentator.indent (report);
		report << "Output vector: ";
		VD.write (report, w);
		report << endl;

		if (!VD.areEqual (v, w))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	factory.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testIdentityApply");

	return ret;
}

/* Test 2: Application of nilpotent linear map to random vectors
 *
 * Generates an index-n nilpotent linear map and applies it to randomly
 * generated vectors n times. Generates the vectors so that they are guaranteed
 * to be outside the matrix's index- n-1 subspace; tests to make sure the
 * vectors do not prematurely reach zero.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 *
 * Return true on success and false on failure
 */

template <class Row, class Field, class Vector>
static bool testNilpotentApply (Field &F, const char *text, VectorFactory<Vector> &factory) 
{
	typedef SparseMatrix0 <Field, Vector, Row> Blackbox;

	ostringstream str;
	str << "Testing nilpotent apply (" << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testNilpotentApply", factory.m ());

	bool ret = true;
	bool even, iter_passed;

	StandardBasisFactory<Field, Row> f1 (F, factory.n ());
	Row tmp;
	f1.next (tmp);  // Small trick: pull the first vector out to shift elements up one row
	Blackbox A (F, f1);

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Matrix:" << endl;
	A.write (report, Blackbox::FORMAT_PRETTY);

	size_t j;
	NonzeroRandIter<Field> r (F, typename Field::RandIter (F));

	VectorDomain<Field> VD (F);
	Vector v, w;

	VectorWrapper::ensureDim (v, factory.n ());
	VectorWrapper::ensureDim (w, factory.n ());

	while (factory) {
		commentator.startIteration (factory.j ());

		iter_passed = true;
		even = false;

		factory.next (v);

		// Make sure last element is nonzero
		r.random (VectorWrapper::ref<Field> (v, factory.n () - 1));

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		VD.write (report, v);
		report << endl;

		commentator.start ("Applying vectors");

		for (j = 0; j < factory.n () - 1; j++, even = !even)
			if (even)
				A.apply (v, w);
			else
				A.apply (w, v);

		commentator.stop ("Done");

		commentator.indent (report);
		report << "A^(n-1) v:     ";
		VD.write (report, even ? w : v);
		report << endl;

		if (VD.isZero (even ? w : v)) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: A^(n-1) v is prematurely zero" << endl;
		}

		if (even)
			A.apply (v, w);
		else
			A.apply (w, v);

		commentator.indent (report);
		report << "A^n v:         ";
		VD.write (report, even ? v : w);
		report << endl;

		if (!VD.isZero (even ? v : w))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: A^n v is non-zero" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	factory.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testNilpotentApply");

	return ret;
}

/* Test 3: Random apply to sparse matrix of K nonzero elements per row
 *
 * Generates a random sparse matrix with K nonzero elements per row and applies
 * it to the vectors {e_i | i=1..n} to test whether the output matches the ith
 * column of the matrix.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * K - Number of nonzero elements per row
 *
 * Return true on success and false on failure
 */

template <class Vector, class Row, class Field>
bool testRandomApply1 (Field &F, const char *text, unsigned int iterations, VectorFactory<Row> &A_factory) 
{
	typedef SparseMatrix0 <Field, Vector, Row> Blackbox;

	ostringstream str;
	str << "Testing sparse random apply (1, " << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testRandomApply1", iterations);

	bool ret = true;
	bool iter_passed;

	size_t i, k;

	VectorDomain<Field> VD (F);

	StandardBasisFactory<Field, Vector> factory (F, A_factory.n ());
	Vector e_j, w;

	VectorWrapper::ensureDim (e_j, A_factory.n ());
	VectorWrapper::ensureDim (w, A_factory.m ());

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		iter_passed = true;

		Blackbox A (F, A_factory);
		A_factory.reset ();

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.write (report, Blackbox::FORMAT_PRETTY);

		factory.reset ();

		while (factory) {
			factory.next (e_j);

			A.apply (w, e_j);

			for (k = 0; k < A_factory.m (); k++)
				if (!F.areEqual (A.getEntry (k, factory.j () - 1), VectorWrapper::constRef<Field> (w, k)))
					ret = iter_passed = false;

			commentator.indent (report);
			report << "Output vector " << factory.j () << ": ";
			VD.write (report, w) << endl;
		}

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Output vectors were incorrect" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomApply1");

	return ret;
}

/* Test 5: Random apply to sparse matrix of K nonzero elements per row
 *
 * Generates a random sparse matrix with K nonzero elements per row and applies
 * it to the vectors (1,1,...,1) to test whether the output matches the sum of
 * the input's columns
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * K - Number of nonzero elements per row
 *
 * Return true on success and false on failure
 */

template <class Vector, class Row, class Field>
bool testRandomApply2 (Field &F, const char *text, unsigned int iterations, VectorFactory<Row> &A_factory) 
{
	typedef SparseMatrix0 <Field, Vector, Row> Blackbox;

	ostringstream str;
	str << "Testing sparse random apply (2, " << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testRandomApply2", iterations);

	bool ret = true;
	bool iter_passed;

	size_t i, j, k;

	VectorDomain<Field> VD (F);
	typename Field::RandIter r (F);
	typename Field::Element sum;

	integer c;
	long width;

	F.characteristic (c);
	width = logp (c, 10) + 1;

	Vector v, w;

	VectorWrapper::ensureDim (v, A_factory.n ());
	VectorWrapper::ensureDim (w, A_factory.m ());

	for (k = 0; k < A_factory.n (); k++)
		F.init (VectorWrapper::ref<Field> (v, k), 1);

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		iter_passed = true;

		Blackbox A (F, A_factory);

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.write (report, Blackbox::FORMAT_PRETTY);

		A.apply (w, v);

		for (j = 0; j < A_factory.m (); j++) {
			F.init (sum, 0);

			for (k = 0; k < A_factory.n (); k++)
				F.addin (sum, A.getEntry (j, k));

			if (!F.areEqual (sum, VectorWrapper::constRef<Field> (w, j)))
				ret = iter_passed = false;
		}

		commentator.indent (report);
		report << "Output vector: ";
		VD.write (report, w) << endl;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Output vector was incorrect" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomApply2");

	return ret;
}

/* Test 6: Random transpose
 *
 * Construct a random sparse matrix and check that the output of applyTranspose
 * is consistent with the output of apply
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random vectors to which to apply matrix
 *
 * Return true on success and false on failure
 */

template <class Row, class Field, class Vector>
static bool testRandomTranspose (Field                 &F,
				 const char            *text,
				 VectorFactory<Row>    &A_factory,
				 VectorFactory<Vector> &factory1,
				 VectorFactory<Vector> &factory2) 
{
	typedef SparseMatrix0 <Field, Vector, Row> Blackbox;

	ostringstream str;
	str << "Testing random transpose (" << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testRandomTranspose", factory1.m ());

	Blackbox A (F, A_factory);
	A_factory.reset ();

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix:" << endl;
	A.write (report, Blackbox::FORMAT_PRETTY);

	bool ret = testTranspose (F, A, factory1, factory2);

	factory1.reset ();
	factory2.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomTranspose");

	return ret;
}

/* Test 7: Linearity
 *
 * Compute a random sparse matrix and check its linearity
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random vectors to which to apply matrix
 *
 * Return true on success and false on failure
 */

template <class Row, class Field, class Vector>
static bool testRandomLinearity (Field                 &F,
				 const char            *text,
				 VectorFactory<Row>    &A_factory,
				 VectorFactory<Vector> &factory1,
				 VectorFactory<Vector> &factory2) 
{
	typedef SparseMatrix0 <Field, Vector, Row> Blackbox;

	ostringstream str;
	str << "Testing linearity (" << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testRandomLinearity", factory1.m ());

	Blackbox A (F, A_factory);
	A_factory.reset ();

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix:" << endl;
	A.write (report, Blackbox::FORMAT_PRETTY);

	bool ret = testLinearity (F, A, factory1, factory2);

	factory1.reset ();
	factory2.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomLinearity");

	return ret;
}

template <class Field, class Vector, class Row>
bool runSparseMatrixTestsByVector (const Field           &F,
				   const char            *desc,
				   unsigned int           iterations,
				   VectorFactory<Vector> &v_factory1,
				   VectorFactory<Vector> &v_factory2,
				   VectorFactory<Row>    &A_factory) 
{
	bool pass = true;

	ostringstream str;

	str << "Testing " << desc << " sparse matrix" << ends;
	commentator.start (str.str ().c_str (), "runSparseMatrixTestsByVector", 6);

	if (!testIdentityApply<Row>   (F, desc, v_factory1))                        pass = false; commentator.progress ();
	if (!testNilpotentApply<Row>  (F, desc, v_factory1))                        pass = false; commentator.progress ();
	if (!testRandomApply1<Vector> (F, desc, iterations, A_factory))             pass = false; commentator.progress ();
	if (!testRandomApply2<Vector> (F, desc, iterations, A_factory))             pass = false; commentator.progress ();
	if (!testRandomTranspose      (F, desc, A_factory, v_factory1, v_factory2)) pass = false; commentator.progress ();
	if (!testRandomLinearity      (F, desc, A_factory, v_factory1, v_factory2)) pass = false; commentator.progress ();

	commentator.stop (MSG_STATUS (pass), (const char *) 0, "runSparseMatrixTests");

	return pass;
}

template <class Field, class Row>
bool runSparseMatrixTests (const Field        &F,
			   const char         *desc,
			   int                 iterations,
			   VectorFactory<Row> &A_factory) 
{
	bool pass = true;

	ostringstream str1, str2, str3, str4, str5;

	str1 << "Testing sparse matrix with row type " << desc << ends;
	commentator.start (str1.str ().c_str (), "runSparseMatrixTests", 4);

	str2 << desc << "/dense" << ends;
	str3 << desc << "/sparse sequence" << ends;
	str4 << desc << "/sparse associative" << ends;
	str5 << desc << "/sparse parallel" << ends;

	RandomDenseVectorFactory<Field>     dense_factory1 (F, A_factory.n (), iterations);
	RandomDenseVectorFactory<Field>     dense_factory2 (F, A_factory.m (), iterations);
	RandomSparseSeqVectorFactory<Field> sparse_seq_factory1 (F, A_factory.n (), A_factory.n () / 10, iterations);
	RandomSparseSeqVectorFactory<Field> sparse_seq_factory2 (F, A_factory.m (), A_factory.m () / 10, iterations);
	RandomSparseMapVectorFactory<Field> sparse_map_factory1 (F, A_factory.n (), A_factory.n () / 10, iterations);
	RandomSparseMapVectorFactory<Field> sparse_map_factory2 (F, A_factory.m (), A_factory.m () / 10, iterations);
	RandomSparseParVectorFactory<Field> sparse_par_factory1 (F, A_factory.n (), A_factory.n () / 10, iterations);
	RandomSparseParVectorFactory<Field> sparse_par_factory2 (F, A_factory.m (), A_factory.m () / 10, iterations);

	if (!runSparseMatrixTestsByVector (F, str2.str ().c_str (), iterations,
					   dense_factory1, dense_factory2, A_factory))
		pass = false;
	commentator.progress ();
	if (!runSparseMatrixTestsByVector (F, str2.str ().c_str (), iterations,
					   sparse_seq_factory1, sparse_seq_factory2, A_factory))
		pass = false;
	commentator.progress ();
	if (!runSparseMatrixTestsByVector (F, str2.str ().c_str (), iterations,
					   sparse_map_factory1, sparse_map_factory2, A_factory))
		pass = false;
	commentator.progress ();
	if (!runSparseMatrixTestsByVector (F, str2.str ().c_str (), iterations,
					   sparse_par_factory1, sparse_par_factory2, A_factory))
		pass = false;
	commentator.progress ();

	commentator.stop (MSG_STATUS (pass), (const char *) 0, "runSparseMatrixTests");

	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static size_t m = 10;
	static integer q = 101;
	static int iterations = 100;
	static int k = 3;
	static int N = 20;

	static Argument args[] = {
		{ 'n', "-n N", "Set column dimension of test matrices to N (default 10)",            TYPE_INT,     &n },
		{ 'm', "-m M", "Set row dimension of test matrices to M (default 10)",               TYPE_INT,     &m },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 101)",                 TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",                   TYPE_INT,     &iterations },
		{ 'k', "-k K", "K nonzero Elements per row in sparse random apply test (default 3)", TYPE_INT,     &k },
		{ 'N', "-N N", "N nonzero Elements in sparse random apply test (default 20)",        TYPE_INT,     &N }
	};

	typedef	Modular<uint32> Field;
	typedef Field::Element  Element;

	typedef std::vector <pair <size_t, Element> > SeqRow;
	typedef std::map <size_t, Element> MapRow;
	typedef std::pair <std::vector<size_t>, std::vector<Element> > ParRow;

	typedef std::vector <Element> DenseVector;
	typedef std::vector <pair <size_t, Element> > SparseSeqVector;
	typedef std::map <size_t, Element> SparseMapVector;
	typedef std::pair <std::vector<size_t>, std::vector<Element> > SparseParVector;

	parseArguments (argc, argv, args);
	Field F (q);

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	cout << "Sparse matrix black box test suite" << endl << endl;

	NonzeroRandIter<Field> r (F, Field::RandIter (F));

	RandomSparseSeqVectorFactory<Field, NonzeroRandIter<Field> > factory1 (F, r, n, k, m);
	RandomSparseMapVectorFactory<Field, NonzeroRandIter<Field> > factory2 (F, r, n, k, m);
	RandomSparseParVectorFactory<Field, NonzeroRandIter<Field> > factory3 (F, r, n, k, m);

	if (!runSparseMatrixTests (F, "sparse sequence",    iterations, factory1)) pass = false;
	if (!runSparseMatrixTests (F, "sparse associative", iterations, factory2)) pass = false;
	if (!runSparseMatrixTests (F, "sparse parallel",    iterations, factory3)) pass = false;

	return pass ? 0 : -1;
}
