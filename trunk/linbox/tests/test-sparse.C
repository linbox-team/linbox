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

template <class Row, class Vector, class Field>
bool testRandomApply1 (Field &F, const char *text, size_t n, size_t iterations, size_t K) 
{
	typedef SparseMatrix0 <Field, Vector, Row> Blackbox;

	ostringstream str;
	str << "Testing sparse random apply (1, " << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testRandomApply1", iterations);

	bool ret = true;
	bool iter_passed;

	size_t i, j, k;

	typename Field::RandIter r (F);
	VectorDomain<Field> VD (F);

	StandardBasisFactory<Field, Vector> factory (F, n);
	Vector e_j, w;

	VectorWrapper::ensureDim (e_j, n);
	VectorWrapper::ensureDim (w, n);

	if (K > n) K = n;

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		iter_passed = true;

		Blackbox A (F, n, n);

		for (j = 0; j < n; j++) {
			for (k = 0; k < K; k++) {
				size_t l;

				do
					l = rand () % n;
				while (!F.isZero (A.getEntry (j, l)));

				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
					<< __FUNCTION__ << ": now ready for non-const version" << endl;

				r.random (A.refEntry (j, l));
			}
		}

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.write (report, Blackbox::FORMAT_PRETTY);

		factory.reset ();

		while (factory) {
			factory.next (e_j);

			A.apply (w, e_j);

			for (k = 0; k < n; k++)
				if (!F.areEqual (A.getEntry (k, factory.j () - 1), VectorWrapper::constRef<Field> (w, k)))
					ret = iter_passed = false;

			commentator.indent (report);
			report << "Output vector " << factory.j () << ": ";
			VD.write (report, w);
			report << endl;
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

/* Test 4: Random apply to sparse matrix of N nonzero elements
 *
 * Generates a random sparse matrix with N nonzero elements and applies it to
 * the vectors {e_i | i=1..n} to test whether the output matches the ith column
 * of the matrix.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * N - Number of nonzero elements
 *
 * Return true on success and false on failure
 */

template <class Row, class Vector, class Field>
bool testRandomApply2 (Field &F, const char *text, size_t n, size_t iterations, size_t N) 
{
	typedef SparseMatrix0 <Field, Vector, Row> Blackbox;

	ostringstream str;
	str << "Testing sparse random apply (2, " << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testRandomApply2", iterations);

	bool ret = true;
	bool iter_passed;

	size_t i, k;

	typename Field::RandIter r (F);
	VectorDomain<Field> VD (F);

	StandardBasisFactory<Field, Vector> factory (F, n);
	Vector e_j, w;

	VectorWrapper::ensureDim (e_j, n);
	VectorWrapper::ensureDim (w, n);

	if (N > n * n) N = n * n;

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		iter_passed = true;

		Blackbox A (F, n, n);

		for (k = 0; k < N; k++) {
			size_t l1, l2;

			do {
				l1 = rand () % n;
				l2 = rand () % n;
			} while (!F.isZero (A.getEntry (l1, l2)));

			r.random (A.refEntry (l1, l2));
		}

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.write (report, Blackbox::FORMAT_PRETTY);

		while (factory) {
			factory.next (e_j);

			A.apply (w, e_j);

			for (k = 0; k < n; k++)
				if (!F.areEqual (A.getEntry (k, factory.j () - 1), VectorWrapper::constRef<Field> (w, k)))
					ret = iter_passed = false;

			commentator.indent (report);
			report << "Output vector " << factory.j () << ": ";
			VD.write (report, w);
			report << endl;
		}

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Output vectors were incorrect" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomApply2");

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

template <class Row, class Vector, class Field>
bool testRandomApply3 (Field &F, const char *text, size_t n, size_t iterations, size_t K) 
{
	typedef SparseMatrix0 <Field, Vector, Row> Blackbox;

	ostringstream str;
	str << "Testing sparse random apply (3, " << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testRandomApply3", iterations);

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

	if (K > n) K = n;

	Vector v, w;

	VectorWrapper::ensureDim (v, n);
	VectorWrapper::ensureDim (w, n);

	for (k = 0; k < n; k++)
		F.init (VectorWrapper::ref<Field> (v, k), 1);

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		iter_passed = true;

		Blackbox A (F, n, n);

		for (j = 0; j < n; j++) {
			for (k = 0; k < K; k++) {
				size_t l;

				do
					l = rand () % n;
				while (!F.isZero (A.getEntry (j, l)));

				r.random (A.refEntry (j, l));
			}
		}

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.write (report, Blackbox::FORMAT_PRETTY);

		A.apply (w, v);

		for (j = 0; j < n; j++) {
			F.init (sum, 0);

			for (k = 0; k < n; k++)
				F.addin (sum, A.getEntry (j, k));

			if (!F.areEqual (sum, VectorWrapper::constRef<Field> (w, j)))
				ret = iter_passed = false;
		}

		commentator.indent (report);
		report << "Output vector: ";
		VD.write (report, w);
		report << endl;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Output vector was incorrect" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomApply3");

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
				 size_t                 K,
				 VectorFactory<Vector> &factory1,
				 VectorFactory<Vector> &factory2) 
{
	typedef SparseMatrix0 <Field, Vector, Row> Blackbox;

	ostringstream str;
	str << "Testing random transpose (" << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testRandomTranspose", factory1.m ());

	Blackbox A (F, factory1.n (), factory2.n ());
	size_t j, k;
	typename Field::RandIter r (F);

	for (j = 0; j < factory1.n (); j++) {
		for (k = 0; k < K; k++) {
			size_t l;

			do
				l = rand () % factory2.n ();
			while (!F.isZero (A.getEntry (j, l)));

			r.random (A.refEntry (j, l));
		}
	}

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
				 size_t                 K,
				 VectorFactory<Vector> &factory1,
				 VectorFactory<Vector> &factory2) 
{
	typedef SparseMatrix0 <Field, Vector, Row> Blackbox;

	ostringstream str;
	str << "Testing linearity (" << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testRandomLinearity", factory1.m ());

	Blackbox A (F, factory1.n (), factory1.n ());
	size_t j, k;
	typename Field::RandIter r (F);

	for (j = 0; j < factory1.n (); j++) {
		for (k = 0; k < K; k++) {
			size_t l;

			do
				l = rand () % factory1.n ();
			while (!F.isZero (A.getEntry (j, l)));

			r.random (A.refEntry (j, l));
		}
	}

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	report << "Input matrix:" << endl;
	A.write (report, Blackbox::FORMAT_PRETTY);

	bool ret = testLinearity (F, A, factory1, factory2);

	factory1.reset ();
	factory2.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomLinearity");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 101;
	static int iterations = 100;
	static int k = 3;
	static int N = 20;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)",                 TYPE_INT,     &n },
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

	srand (time (NULL));

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);

	cout << "Sparse matrix black box test suite" << endl << endl;

	RandomDenseVectorFactory<Field> factory1 (F, n, iterations);
	RandomSparseSeqVectorFactory<Field> factory2 (F, n, n / 10, iterations);
	RandomSparseMapVectorFactory<Field> factory3 (F, n, n / 10, iterations);
	RandomSparseParVectorFactory<Field> factory4 (F, n, n / 10, iterations);

	RandomDenseVectorFactory<Field, NonzeroRandIter<Field> > factory5 (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), n, iterations);
	RandomSparseSeqVectorFactory<Field, NonzeroRandIter<Field> > factory6 (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), n, n / 10, iterations);
	RandomSparseMapVectorFactory<Field, NonzeroRandIter<Field> > factory7 (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), n, n / 10, iterations);
	RandomSparseParVectorFactory<Field, NonzeroRandIter<Field> > factory8 (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), n, n / 10, iterations);

	if (!testIdentityApply<SeqRow>   (F, "sparse sequence/dense",                 factory1)) pass = false;
	if (!testIdentityApply<SeqRow>   (F, "sparse sequence/sparse sequence",       factory2)) pass = false;
	if (!testIdentityApply<SeqRow>   (F, "sparse sequence/sparse associative",    factory3)) pass = false;
	if (!testIdentityApply<SeqRow>   (F, "sparse sequence/sparse parallel",       factory4)) pass = false;
	if (!testIdentityApply<MapRow>   (F, "sparse associative/dense",              factory1)) pass = false;
	if (!testIdentityApply<MapRow>   (F, "sparse associative/sparse sequence",    factory2)) pass = false;
	if (!testIdentityApply<MapRow>   (F, "sparse associative/sparse associative", factory3)) pass = false;
	if (!testIdentityApply<MapRow>   (F, "sparse associative/sparse parallel",    factory4)) pass = false;
	if (!testIdentityApply<ParRow>   (F, "sparse parallel/dense",                 factory1)) pass = false;
	if (!testIdentityApply<ParRow>   (F, "sparse parallel/sparse sequence",       factory2)) pass = false;
	if (!testIdentityApply<ParRow>   (F, "sparse parallel/sparse associative",    factory3)) pass = false;
	if (!testIdentityApply<ParRow>   (F, "sparse parallel/sparse parallel",       factory4)) pass = false;
	if (!testNilpotentApply<SeqRow>  (F, "sparse sequence/dense",                 factory5)) pass = false;
	if (!testNilpotentApply<SeqRow>  (F, "sparse sequence/sparse sequence",       factory6)) pass = false;
	if (!testNilpotentApply<SeqRow>  (F, "sparse sequence/sparse associatve",     factory7)) pass = false;
	if (!testNilpotentApply<SeqRow>  (F, "sparse sequence/sparse parallel",       factory8)) pass = false;
	if (!testNilpotentApply<MapRow>  (F, "sparse associative/dense",              factory5)) pass = false;
	if (!testNilpotentApply<MapRow>  (F, "sparse associative/sparse sequence",    factory6)) pass = false;
	if (!testNilpotentApply<MapRow>  (F, "sparse associative/sparse associative", factory7)) pass = false;
	if (!testNilpotentApply<MapRow>  (F, "sparse associative/sparse parallel",    factory8)) pass = false;
	if (!testNilpotentApply<ParRow>  (F, "sparse parallel/dense",                 factory5)) pass = false;
	if (!testNilpotentApply<ParRow>  (F, "sparse parallel/sparse sequence",       factory6)) pass = false;
	if (!testNilpotentApply<ParRow>  (F, "sparse parallel/sparse associative",    factory7)) pass = false;
	if (!testNilpotentApply<ParRow>  (F, "sparse parallel/sparse parallel",       factory8)) pass = false;
	if (!testRandomApply1<SeqRow, DenseVector>     (F, "sparse sequence/dense",                 n, iterations, k)) pass = false;
	if (!testRandomApply1<SeqRow, SparseSeqVector> (F, "sparse sequence/sparse sequence",       n, iterations, k)) pass = false;
	if (!testRandomApply1<SeqRow, SparseMapVector> (F, "sparse sequence/sparse associative",    n, iterations, k)) pass = false;
	if (!testRandomApply1<SeqRow, SparseParVector> (F, "sparse sequence/sparse parallel",       n, iterations, k)) pass = false;
	if (!testRandomApply1<MapRow, DenseVector>     (F, "sparse associative/dense",              n, iterations, k)) pass = false;
	if (!testRandomApply1<MapRow, SparseSeqVector> (F, "sparse associative/sparse sequence",    n, iterations, k)) pass = false;
	if (!testRandomApply1<MapRow, SparseMapVector> (F, "sparse associative/sparse associative", n, iterations, k)) pass = false;
	if (!testRandomApply1<MapRow, SparseParVector> (F, "sparse associative/sparse parallel",    n, iterations, k)) pass = false;
	if (!testRandomApply1<ParRow, DenseVector>     (F, "sparse parallel/dense",                 n, iterations, k)) pass = false;
	if (!testRandomApply1<ParRow, SparseSeqVector> (F, "sparse parallel/sparse sequence",       n, iterations, k)) pass = false;
	if (!testRandomApply1<ParRow, SparseMapVector> (F, "sparse parallel/sparse associative",    n, iterations, k)) pass = false;
	if (!testRandomApply1<ParRow, SparseParVector> (F, "sparse parallel/sparse parallel",       n, iterations, k)) pass = false;
	if (!testRandomApply2<SeqRow, DenseVector>     (F, "sparse sequence/dense",                 n, iterations, N)) pass = false;
	if (!testRandomApply2<SeqRow, SparseSeqVector> (F, "sparse sequence/sparse sequence",       n, iterations, N)) pass = false;
	if (!testRandomApply2<SeqRow, SparseMapVector> (F, "sparse sequence/sparse associative",    n, iterations, N)) pass = false;
	if (!testRandomApply2<SeqRow, SparseParVector> (F, "sparse sequence/sparse parallel",       n, iterations, N)) pass = false;
	if (!testRandomApply2<MapRow, DenseVector>     (F, "sparse associative/dense",              n, iterations, N)) pass = false;
	if (!testRandomApply2<MapRow, SparseSeqVector> (F, "sparse associative/sparse sequence",    n, iterations, N)) pass = false;
	if (!testRandomApply2<MapRow, SparseMapVector> (F, "sparse associative/sparse associative", n, iterations, N)) pass = false;
	if (!testRandomApply2<MapRow, SparseParVector> (F, "sparse associative/sparse parallel",    n, iterations, N)) pass = false;
	if (!testRandomApply2<ParRow, DenseVector>     (F, "sparse parallel/dense",                 n, iterations, N)) pass = false;
	if (!testRandomApply2<ParRow, SparseSeqVector> (F, "sparse parallel/sparse sequence",       n, iterations, N)) pass = false;
	if (!testRandomApply2<ParRow, SparseMapVector> (F, "sparse parallel/sparse associative",    n, iterations, N)) pass = false;
	if (!testRandomApply2<ParRow, SparseParVector> (F, "sparse parallel/sparse parallel",       n, iterations, N)) pass = false;
	if (!testRandomApply3<SeqRow, DenseVector>     (F, "sparse sequence/dense",                 n, iterations, k)) pass = false;
	if (!testRandomApply3<SeqRow, SparseSeqVector> (F, "sparse sequence/sparse sequence",       n, iterations, k)) pass = false;
	if (!testRandomApply3<SeqRow, SparseMapVector> (F, "sparse sequence/sparse associative",    n, iterations, k)) pass = false;
	if (!testRandomApply3<SeqRow, SparseParVector> (F, "sparse sequence/sparse parallel",       n, iterations, k)) pass = false;
	if (!testRandomApply3<MapRow, DenseVector>     (F, "sparse associative/dense",              n, iterations, k)) pass = false;
	if (!testRandomApply3<MapRow, SparseSeqVector> (F, "sparse associative/sparse sequence",    n, iterations, k)) pass = false;
	if (!testRandomApply3<MapRow, SparseMapVector> (F, "sparse associative/sparse associative", n, iterations, k)) pass = false;
	if (!testRandomApply3<MapRow, SparseParVector> (F, "sparse associative/sparse parallel",    n, iterations, k)) pass = false;
	if (!testRandomApply3<ParRow, DenseVector>     (F, "sparse parallel/dense",                 n, iterations, k)) pass = false;
	if (!testRandomApply3<ParRow, SparseSeqVector> (F, "sparse parallel/sparse sequence",       n, iterations, k)) pass = false;
	if (!testRandomApply3<ParRow, SparseMapVector> (F, "sparse parallel/sparse associative",    n, iterations, k)) pass = false;
	if (!testRandomApply3<ParRow, SparseParVector> (F, "sparse parallel/sparse parallel",       n, iterations, k)) pass = false;
	if (!testRandomTranspose<SeqRow> (F, "sparse sequence/dense",                 n, factory1, factory5)) pass = false;
	if (!testRandomTranspose<SeqRow> (F, "sparse sequence/sparse sequence",       n, factory2, factory6)) pass = false;
	if (!testRandomTranspose<SeqRow> (F, "sparse sequence/sparse associative",    n, factory3, factory7)) pass = false;
	if (!testRandomTranspose<SeqRow> (F, "sparse sequence/sparse parallel",       n, factory4, factory8)) pass = false;
	if (!testRandomTranspose<MapRow> (F, "sparse associative/dense",              n, factory1, factory5)) pass = false;
	if (!testRandomTranspose<MapRow> (F, "sparse associative/sparse sequence",    n, factory2, factory6)) pass = false;
	if (!testRandomTranspose<MapRow> (F, "sparse associative/sparse associative", n, factory3, factory7)) pass = false;
	if (!testRandomTranspose<MapRow> (F, "sparse associative/sparse parallel",    n, factory4, factory8)) pass = false;
	if (!testRandomTranspose<ParRow> (F, "sparse parallel/dense",                 n, factory1, factory5)) pass = false;
	if (!testRandomTranspose<ParRow> (F, "sparse parallel/sparse sequence",       n, factory2, factory6)) pass = false;
	if (!testRandomTranspose<ParRow> (F, "sparse parallel/sparse associative",    n, factory3, factory7)) pass = false;
	if (!testRandomTranspose<ParRow> (F, "sparse parallel/sparse parallel",       n, factory4, factory8)) pass = false;
	if (!testRandomLinearity<SeqRow> (F, "sparse sequence/dense",                 k, factory1, factory5)) pass = false;
	if (!testRandomLinearity<SeqRow> (F, "sparse sequence/sparse sequence",       k, factory2, factory6)) pass = false;
	if (!testRandomLinearity<SeqRow> (F, "sparse sequence/sparse associative",    k, factory3, factory7)) pass = false;
	if (!testRandomLinearity<SeqRow> (F, "sparse sequence/sparse parallel",       k, factory4, factory8)) pass = false;
	if (!testRandomLinearity<MapRow> (F, "sparse associative/dense",              k, factory1, factory5)) pass = false;
	if (!testRandomLinearity<MapRow> (F, "sparse associative/sparse sequence",    k, factory2, factory6)) pass = false;
	if (!testRandomLinearity<MapRow> (F, "sparse associative/sparse associative", k, factory3, factory7)) pass = false;
	if (!testRandomLinearity<MapRow> (F, "sparse associative/sparse parallel",    k, factory4, factory8)) pass = false;
	if (!testRandomLinearity<ParRow> (F, "sparse parallel/dense",                 k, factory1, factory5)) pass = false;
	if (!testRandomLinearity<ParRow> (F, "sparse parallel/sparse sequence",       k, factory2, factory6)) pass = false;
	if (!testRandomLinearity<ParRow> (F, "sparse parallel/sparse associative",    k, factory3, factory7)) pass = false;
	if (!testRandomLinearity<ParRow> (F, "sparse parallel/sparse parallel",       k, factory4, factory8)) pass = false;

	return pass ? 0 : -1;
}
