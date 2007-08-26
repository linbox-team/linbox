/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-trace.C
 * Copyright (C) 2002 Bradford Hovinen
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
#include <vector>
#include <cstdio>

#include "test-common.h"

#include "linbox/util/commentator.h"
#include "linbox/field/modular-int.h"
#include "linbox/solutions/trace.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/dense.h"
#include "linbox/vector/stream.h"

using namespace LinBox;

/* Test 1: Trace of random diagonal matrix
 *
 * Construct a random diagonal matrix and check that its computed trace is the
 * same as the sum of its entries
 *
 * F - Field over which to perform computations
 * stream - Stream that comprises source of diagonal vectors
 *
 * Return true on success and false on failure
 */

template <class Field>
bool testScalarMatrixTrace (const Field &F, size_t n)
{
	bool ret = true;
	commentator.start ("Testing scalar matrix trace", "", 1);
	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "scalarmatrix trace test (using specialization)" << endl;
    	typename Field::Element s, t, th; 
	F.init(s, 2);
	F.init(th, 2*n);
	ScalarMatrix<Field> B(F, n, s);
	trace(t, B);
	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testScalarMatrixTrace");
	if (!F.areEqual(t, th)) {
	report << "bad scalar matrix trace " << t << ", should be " << th << endl;

		return false; 
	} 
	else return true;
}

template <class Field>
bool testSparseMatrixTrace (const Field &F, size_t n)
{
	commentator.start ("Building sparse matrix", "", 1);
	bool ret = true;
	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
    	typename Field::Element s, t, th; 
	F.init(s, 2);
	size_t m = (n > 10 ? 10 : n);
	F.init(th, 2*m);
	SparseMatrix<Field> B(F, n, n);
	for (size_t i = 0; i <  m; ++i)
		for (size_t j = 0; j < m; ++j) 
			B.setEntry(i,j,s);
	commentator.stop ("", "done");
	commentator.start ("Testing sparse matrix trace", "", 1);
	report << "sparse matrix trace test (using specialization)" << endl;
	trace(t, B);
	if (!F.areEqual(t, th)) {
	report << "bad sparse matrix trace " << t << ", should be " << th << endl;

		ret = false; 
	} 
	else ret = true;
	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testSparseMatrixTrace");
	return ret;
}

template <class Field>
static bool testDenseMatrixTrace (const Field &F, size_t n)
{
	bool ret = true;
    	typename Field::Element s, t, th; 
	F.init(s, 2);
	size_t m = (n > 10 ? 10 : n);
	F.init(th, 2*m);
	DenseMatrix<Field> B(F, n, n);
	for (size_t i = 0; i <  m; ++i)
		for (size_t j = 0; j < n; ++j) 
			B.setEntry(i, j, s);
	commentator.start ("Testing dense matrix trace", "", 1);
	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "dense matrix trace test (using specialization)" << endl;
	trace(t, B);
	if (!F.areEqual(t, th)) {
		report << "bad dense matrix trace " << t << ", should be " << th << endl;

		ret = false; 
	} 
	else ret = true;
	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDenseMatrixTrace");
	return ret;
}

template <class Field>
static bool testDiagonalTrace (const Field &F, VectorStream<vector<typename Field::Element> > &stream) 
{
	typedef vector <typename Field::Element> Vector;
	typedef Diagonal <Field> Blackbox;

	commentator.start ("Testing diagonal trace", "testDiagonalTrace", stream.m ());
	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	VectorDomain<Field> VD (F);

	bool ret = true;
	size_t i;

	Vector d;
	typename Field::Element sigma, res;

	VectorWrapper::ensureDim (d, stream.dim ());

	while (stream) {
		commentator.startIteration (stream.j ());

		stream.next (d);

		report << "Input vector:  ";
		VD.write (report, d);
		report << endl;

		F.init (sigma, 0);
		for (i = 0; i < stream.n (); i++)
			F.addin (sigma, VectorWrapper::constRef<Field, Vector> (d, i));

		report << "True trace: ";
		F.write (report, sigma);
		report << endl;

		Blackbox D (F, d);

		trace (res, D);

		report << "Computed trace: ";
		F.write (report, res);
		report << endl;

		if (!F.areEqual (sigma, res)) {
			ret = false;
			report << "ERROR: Computed trace is incorrect" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalTrace");

	return ret;
}

template <class Field>
static bool testComposeTrace (const Field &F, size_t n, VectorStream<vector<typename Field::Element> > &stream)
{
	bool ret = true;
    	typename Field::Element s; 
	F.init(s, 2);

	SparseMatrix<Field> B(F, n, n);
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < n; ++j) 
			B.setEntry(i, j, s);

	VectorDomain<Field> VD (F);
        

	commentator.start ("Testing composed matrix trace", "", 1);
	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "composed trace test (using specialization)" << endl;

	typedef vector <typename Field::Element> Vector;

	size_t i;

	Vector d;
	typename Field::Element sigma, res, th;

	VectorWrapper::ensureDim (d, stream.dim ());

	while (stream) {
		commentator.startIteration (stream.j ());

		stream.next (d);

		report << "Input diagonal:  ";
		VD.write (report, d);
		report << endl;

                B.write (report << "Input dense " , FORMAT_MAPLE) << endl;

		F.init (sigma, 0);
		F.init (th, 0);
		for (i = 0; i < stream.n (); i++) {
			F.addin (sigma, VectorWrapper::constRef<Field, Vector> (d, i));
                        F.axpyin(th, VectorWrapper::constRef<Field, Vector> (d, i), VectorWrapper::constRef<Field, Vector> (d, i) );
                }
                F.mulin(sigma, s);
                F.mulin(th, s);

		report << "True trace: ";
		F.write (report, sigma);
		report << endl;

		Diagonal<Field> D (F, d);

                Compose< Diagonal<Field>, SparseMatrix<Field> > CDB(&D, &B);

		trace (res, CDB);

		report << "Computed trace: ";
		F.write (report, res);
		report << endl;
                
		if (!F.areEqual (sigma, res)) {
			ret = false;
			report << "ERROR: Computed trace is incorrect" << endl;
		}


		report << "True trace: ";
		F.write (report, th);
		report << endl;

                Compose< Compose< Diagonal<Field>, SparseMatrix<Field> >, Diagonal<Field> > CDBD(&CDB, &D);

		trace (res, CDBD);

		report << "Computed trace: ";
		F.write (report, res);
		report << endl;
                
		if (!F.areEqual (th, res)) {
			ret = false;
			report << "ERROR: Computed trace is incorrect" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testComposedMatrixTrace");
	return ret;
}



int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 256;
	//static size_t n = 10;
	static integer q = 101;
	static int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};

	typedef Modular<int> Field;
	typedef vector<Field::Element> Vector;

	parseArguments (argc, argv, args);
	Field F (q);

	commentator.start("Trace solution test suite", "Trace");
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);

	RandomDenseStream<Field, Vector> stream (F, n, iterations);

	if (!testScalarMatrixTrace (F, n)) pass = false;
	if (!testSparseMatrixTrace (F, n)) pass = false;
	if (!testDenseMatrixTrace (F, n)) pass = false;
	if (!testDiagonalTrace (F, stream)) pass = false;
        stream.reset();
	if (!testComposeTrace (F, n, stream)) pass = false;

	commentator.stop("Trace solution test suite");
	return pass ? 0 : -1;
}
