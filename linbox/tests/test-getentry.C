/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-getentry.C
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
#include "linbox/field/modular-int32.h"
#include "linbox/solutions/getentry.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/dense.h"
#include "linbox/vector/stream.h"

using namespace LinBox;

/* Test 1: getEntry of random diagonal matrix
 *
 * Construct a random diagonal matrix and check that its computed getEntry is the
 * correct i,j element
 *
 * F - Field over which to perform computations
 * stream - Stream that comprises source of diagonal vectors
 *
 * Return true on success and false on failure
 */

template <class Field>
bool testScalarMatrixgetEntry (const Field &F, size_t n)
{
	bool ret = true;
	commentator.start ("Testing scalar matrix getEntry", "", 1);
	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "scalarmatrix getEntry test (using specialization)" << endl;
    	typename Field::Element s, t, r, th; 
	F.init(r, 0);
	F.init(s, 2);
	F.init(th, 2);
	ScalarMatrix<Field> B(F, n, s);
	getEntry(t, B, 0, n-1); F.assign(r,t);
		report << "0xn-1" << t << endl;
  	getEntry(t, B, n-1, 0); F.addin(r,t);
		report << "n-1 x  0" << t << endl;
  	getEntry(t, B, 0, 0); F.addin(r,t);
		report << "0 x 0" << t << endl;
  	getEntry(t, B, n-1, n-1); F.addin(r,t);
		report << "n-1 x n-1" << t << endl;
	if (!F.areEqual(t, th)) {
		report << "bad scalar matrix getEntry " << t << ", should be " << th << endl;
		ret= false; 
	} 
	else ret= true;
	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testScalarMatrixgetEntry");
	return ret;
}

template <class Field>
bool testSparseMatrixgetEntry (const Field &F, size_t n)
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
	commentator.start ("Testing sparse matrix getEntry", "", 1);
	report << "sparse matrix getEntry test (using specialization)" << endl;
        ret = true;
	getEntry(t, B, 0, 0);
	if (!F.areEqual(t, s)) {
	report << "bad sparse matrix getEntry 1,1 " << t << ", should be " << s << endl;

		ret = false; 
	} 
	getEntry(t, B, 0, n-1);
	if (!F.areEqual(t, s)) {
	report << "bad sparse matrix getEntry 1,n" << t << ", should be " << s << endl;

		ret = false; 
	} 
	getEntry(t, B, n-1, 0);
	if (!F.areEqual(t, s)) {
	report << "bad sparse matrix getEntry n,1" << t << ", should be " << s << endl;

		ret = false; 
	} 
	getEntry(t, B, n-1, n-1);
	if (!F.areEqual(t, s)) {
	report << "bad sparse matrix getEntry n,n" << t << ", should be " << s << endl;

		ret = false; 
	} 
	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testSparseMatrixgetEntry");
	return ret;
}

template <class Field>
static bool testDenseMatrixgetEntry (const Field &F, size_t n)
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
	commentator.start ("Testing dense matrix getEntry", "", 1);
	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "dense matrix getEntry test (using specialization)" << endl;
        ret = true;
        getEntry(t, B, 0, 0);
	if (!F.areEqual(t, s)) {
		report << "bad dense matrix getEntry 1,1 " << t << ", should be " << s << endl;

		ret = false; 
	} 
        getEntry(t, B, 0, n-1);
	if (!F.areEqual(t, s)) {
		report << "bad dense matrix getEntry 1,n " << t << ", should be " << s << endl;

		ret = false; 
	} 
        getEntry(t, B, n-1, 0);
	if (!F.areEqual(t, s)) {
		report << "bad dense matrix getEntry n,1 " << t << ", should be " << s << endl;

		ret = false; 
	} 
        getEntry(t, B, n-1, n-1);
	if (!F.areEqual(t, s)) {
		report << "bad dense matrix getEntry n,n " << t << ", should be " << s << endl;

		ret = false; 
	} 
	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDenseMatrixgetEntry");
	return ret;
}

template <class Field>
static bool testDiagonalgetEntry (const Field &F, VectorStream<vector<typename Field::Element> > &stream) 
{
	typedef vector <typename Field::Element> Vector;
	typedef Diagonal <Field> Blackbox;

	commentator.start ("Testing diagonal getEntry", "testDiagonalgetEntry", stream.m ());
	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	VectorDomain<Field> VD (F);

	bool ret = true;
	size_t i;

	Vector d;
	typename Field::Element sigma, res, ge;

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

		report << "True getEntry: ";
		F.write (report, sigma);
		report << endl;

		Blackbox D (F, d);


		F.init (res, 0);
		for (i = 0; i < stream.n (); i++)
			F.addin (res, getEntry (ge, D, i, i));
                F.addin(res, getEntry (ge, D, 0,stream.n ()-1));
                F.addin(res, getEntry (ge, D, stream.n ()-1, 0));
               

		report << "Computed getEntry: ";
		F.write (report, res);
		report << endl;

		if (!F.areEqual (sigma, res)) {
			ret = false;
			report << "ERROR: Computed getEntry is incorrect" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalgetEntry");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	//static size_t n = 256;
	static size_t n = 10;
	static integer q = 101;
	static int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 256)", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 101)",  TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",     TYPE_INT,     &iterations },
	};

	typedef Modular<int32> Field;
	typedef vector<Field::Element> Vector;

	parseArguments (argc, argv, args);
	Field F (q);

	cout << endl << "Black box getEntry test suite" << endl;
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);

	RandomDenseStream<Field, Vector> stream (F, n, iterations);

	if (!testScalarMatrixgetEntry (F, n)) pass = false;
	if (!testSparseMatrixgetEntry (F, n)) pass = false;
	if (!testDenseMatrixgetEntry (F, n)) pass = false;
	if (!testDiagonalgetEntry (F, stream)) pass = false;

	return pass ? 0 : -1;
}
