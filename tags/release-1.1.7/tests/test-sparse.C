
/* tests/test-sparse.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/vector/stream.h"
#include "linbox/vector/vector-domain.h"
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
static bool testIdentityApply (Field &F, const char *text, VectorStream<Vector> &stream) 
{
	typedef SparseMatrix <Field, Row> Blackbox;

	ostringstream str;
	str << "Testing identity apply (" << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testIdentityApply", stream.m ());

	bool ret = true;
	bool iter_passed = true;

	VectorDomain<Field> VD (F);
	StandardBasisStream<Field, Row> f1 (F, stream.n ());
	Blackbox A (F, f1);

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Matrix:" << endl;
	A.write (report, FORMAT_PRETTY);

	Vector v, w;

	VectorWrapper::ensureDim (v, stream.n ());
	VectorWrapper::ensureDim (w, stream.n ());

	while (stream) {
		commentator.startIteration (stream.j ());

		iter_passed = true;

		stream.next (v);

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		VD.write (report, v);
		report << endl;

		A.apply (w, v);

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

	stream.reset ();

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
static bool testNilpotentApply (Field &F, const char *text, VectorStream<Vector> &stream) 
{
	typedef SparseMatrix <Field, Row> Blackbox;

	ostringstream str;
	str << "Testing nilpotent apply (" << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testNilpotentApply", stream.m ());

	bool ret = true;
	bool even, iter_passed;

	StandardBasisStream<Field, Row> f1 (F, stream.n ());
	Row tmp;
	f1.next (tmp);  // Small trick: pull the first vector out to shift elements up one row
	Blackbox A (F, f1);

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Matrix:" << endl;
	A.write (report, FORMAT_PRETTY);

	size_t j;
	NonzeroRandIter<Field> r (F, typename Field::RandIter (F));

	VectorDomain<Field> VD (F);
	Vector v, w;

	VectorWrapper::ensureDim (v, stream.n ());
	VectorWrapper::ensureDim (w, stream.n ());

	while (stream) {
		commentator.startIteration (stream.j ());

		iter_passed = true;
		even = false;

		stream.next (v);

		// Make sure last element is nonzero
		r.random (VectorWrapper::ref<Field> (v, stream.n () - 1));

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		VD.write (report, v);
		report << endl;

		commentator.start ("Applying vectors");

		for (j = 0; j < stream.n () - 1; j++, even = !even)
			if (even)
				A.apply (v, w);
			else
				A.apply (w, v);

		commentator.stop ("Done");

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

	stream.reset ();

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
bool testRandomApply1 (Field &F, const char *text, unsigned int iterations, VectorStream<Row> &A_stream) 
{
	typedef SparseMatrix <Field, Row> Blackbox;

	ostringstream str;
	str << "Testing sparse random apply (1, " << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testRandomApply1", iterations);

	bool ret = true;
	bool iter_passed;

	size_t i, k;

	VectorDomain<Field> VD (F);

	StandardBasisStream<Field, Vector> stream (F, A_stream.n ());
	Vector e_j, w;

	VectorWrapper::ensureDim (e_j, A_stream.n ());
	VectorWrapper::ensureDim (w, A_stream.m ());

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		iter_passed = true;

		Blackbox A (F, A_stream);
		A_stream.reset ();

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.write (report, FORMAT_PRETTY);

		stream.reset ();

		while (stream) {
			stream.next (e_j);

			A.apply (w, e_j);

			for (k = 0; k < A_stream.m (); k++)
				if (!F.areEqual (A.getEntry (k, stream.j () - 1), VectorWrapper::constRef<Field> (w, k)))
					ret = iter_passed = false;

			report << "Output vector " << stream.j () << ": ";
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
bool testRandomApply2 (Field &F, const char *text, unsigned int iterations, VectorStream<Row> &A_stream) 
{
	typedef SparseMatrix <Field, Row> Blackbox;

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

	VectorWrapper::ensureDim (v, A_stream.n ());
	VectorWrapper::ensureDim (w, A_stream.m ());

	for (k = 0; k < A_stream.n (); k++)
		F.init (VectorWrapper::ref<Field> (v, k), 1);

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		iter_passed = true;

		Blackbox A (F, A_stream);
		A_stream.reset ();

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.write (report, FORMAT_PRETTY);

		A.apply (w, v);

		for (j = 0; j < A_stream.m (); j++) {
			F.init (sum, 0);

			for (k = 0; k < A_stream.n (); k++)
				F.addin (sum, A.getEntry (j, k));

			if (!F.areEqual (sum, VectorWrapper::constRef<Field> (w, j)))
				ret = iter_passed = false;
		}

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
static bool testRandomTranspose (Field                &F,
				 const char           *text,
				 VectorStream<Row>    &A_stream,
				 VectorStream<Vector> &stream1,
				 VectorStream<Vector> &stream2) 
{
	typedef SparseMatrix <Field, Row> Blackbox;

	ostringstream str;
	str << "Testing random transpose (" << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testRandomTranspose", stream1.m ());

	Blackbox A (F, A_stream);
	A_stream.reset ();

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix:" << endl;
	A.write (report, FORMAT_PRETTY);

	bool ret = testTranspose (F, A, stream1, stream2);

	stream1.reset ();
	stream2.reset ();

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
				 VectorStream<Row>    &A_stream,
				 VectorStream<Vector> &stream1,
				 VectorStream<Vector> &stream2) 
{
	typedef SparseMatrix <Field, Row> Blackbox;

	ostringstream str;
	str << "Testing linearity (" << text << ")" << ends;
	commentator.start (str.str ().c_str (), "testRandomLinearity", stream1.m ());

	Blackbox A (F, A_stream);
	A_stream.reset ();

	ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Input matrix:" << endl;
	A.write (report, FORMAT_PRETTY);

	bool ret = testLinearity (F, A, stream1, stream2);

	stream1.reset ();
	stream2.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomLinearity");

	return ret;
}

template <class Field, class Vector, class Row>
bool runSparseMatrixTestsByVector (const Field           &F,
				   const char            *desc,
				   unsigned int           iterations,
				   VectorStream<Vector> &v_stream1,
				   VectorStream<Vector> &v_stream2,
				   VectorStream<Row>    &A_stream) 
{
	bool pass = true;

	ostringstream str;

	str << "Testing " << desc << " sparse matrix" << ends;
	commentator.start (str.str ().c_str (), "runSparseMatrixTestsByVector", 6);

	if (!testIdentityApply<Row>   (F, desc, v_stream1))                        pass = false; commentator.progress ();
	if (!testNilpotentApply<Row>  (F, desc, v_stream1))                        pass = false; commentator.progress ();
	if (!testRandomApply1<Vector> (F, desc, iterations, A_stream))             pass = false; commentator.progress ();
	if (!testRandomApply2<Vector> (F, desc, iterations, A_stream))             pass = false; commentator.progress ();
	if (!testRandomTranspose      (F, desc, A_stream, v_stream1, v_stream2)) pass = false; commentator.progress ();
	if (!testRandomLinearity      (F, desc, A_stream, v_stream1, v_stream2)) pass = false; commentator.progress ();

	commentator.stop (MSG_STATUS (pass), (const char *) 0, "runSparseMatrixTests");

	return pass;
}

template <class Field, class Row>
bool runSparseMatrixTests (const Field       &F,
			   const char        *desc,
			   int                iterations,
			   VectorStream<Row> &A_stream) 
{
	typedef std::vector <typename Field::Element> DenseVector;
	typedef std::vector <pair <size_t, typename Field::Element> > SparseSeqVector;
	typedef std::map <size_t, typename Field::Element> SparseMapVector;
	typedef std::pair <std::vector<size_t>, std::vector<typename Field::Element> > SparseParVector;

	bool pass = true;

	ostringstream str1, str2, str3, str4, str5;

	str1 << "Testing sparse matrix with " << desc << " row type" << ends;
	commentator.start (str1.str ().c_str (), "runSparseMatrixTests", 4);

	str2 << desc << "/dense" << ends;
	str3 << desc << "/sparse sequence" << ends;
	str4 << desc << "/sparse associative" << ends;
	str5 << desc << "/sparse parallel" << ends;

	RandomDenseStream<Field, DenseVector>     dense_stream1 (F, A_stream.n (), iterations);
	RandomDenseStream<Field, DenseVector>     dense_stream2 (F, A_stream.m (), iterations);
#if 0
	RandomSparseStream<Field, SparseSeqVector> sparse_seq_stream1 (F, 0.1, A_stream.n (), iterations);
	RandomSparseStream<Field, SparseSeqVector> sparse_seq_stream2 (F, 0.1, A_stream.m (), iterations);
	RandomSparseStream<Field, SparseMapVector> sparse_map_stream1 (F, 0.1, A_stream.n (), iterations);
	RandomSparseStream<Field, SparseMapVector> sparse_map_stream2 (F, 0.1, A_stream.m (), iterations);
	RandomSparseStream<Field, SparseParVector> sparse_par_stream1 (F, 0.1, A_stream.n (), iterations);
	RandomSparseStream<Field, SparseParVector> sparse_par_stream2 (F, 0.1, A_stream.m (), iterations);
#endif

	if (!runSparseMatrixTestsByVector (F, str2.str ().c_str (), iterations,
					   dense_stream1, dense_stream2, A_stream))
		pass = false;
#if 0
	commentator.progress ();
	if (!runSparseMatrixTestsByVector (F, str2.str ().c_str (), iterations,
					   sparse_seq_stream1, sparse_seq_stream2, A_stream))
		pass = false;
	commentator.progress ();
	if (!runSparseMatrixTestsByVector (F, str2.str ().c_str (), iterations,
					   sparse_map_stream1, sparse_map_stream2, A_stream))
		pass = false;
	commentator.progress ();
	if (!runSparseMatrixTestsByVector (F, str2.str ().c_str (), iterations,
					   sparse_par_stream1, sparse_par_stream2, A_stream))
		pass = false;
	commentator.progress ();
#endif

	commentator.stop (MSG_STATUS (pass), (const char *) 0, "runSparseMatrixTests");

	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static size_t m = 10;
	static integer q = 101;
	static int iterations = 1;
	static int k = 3;
	static int N = 20;

	static Argument args[] = {
		{ 'n', "-n N", "Set column dimension of test matrices to N.", TYPE_INT,     &n },
		{ 'm', "-m M", "Set row dimension of test matrices to M.", TYPE_INT,     &m },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 'k', "-k K", "K nonzero Elements per row in sparse random apply test.", TYPE_INT,     &k },
		{ 'N', "-N N", "N nonzero Elements in sparse random apply test.", TYPE_INT,     &N },
		{ '\0' }
	};
	parseArguments (argc, argv, args);

	typedef	Modular<uint32> Field;
	typedef Field::Element  Element;

	typedef std::vector <Element> DenseVector;
	typedef std::vector <pair <size_t, Element> > SparseSeqVector;
	typedef std::map <size_t, Element> SparseMapVector;
	typedef std::pair <std::vector<size_t>, std::vector<Element> > SparseParVector;

	Field F (q);

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	commentator.start("Sparse matrix black box test suite", "Sparse");

	RandomSparseStream<Field, SparseSeqVector>
		stream1 (F, (double) k / (double) n, n, m);
	RandomSparseStream<Field, SparseMapVector>
		stream2 (F, (double) k / (double) n, n, m);
	RandomSparseStream<Field, SparseParVector>
		stream3 (F, (double) k / (double) n, n, m);

	if (!runSparseMatrixTests (F, "sparse sequence",    iterations, stream1)) pass = false;
	if (!runSparseMatrixTests (F, "sparse associative", iterations, stream2)) pass = false;
	if (!runSparseMatrixTests (F, "sparse parallel",    iterations, stream3)) pass = false;

	commentator.stop("Sparse matrix black box test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
