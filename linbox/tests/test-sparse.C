/* -*- mode: c; style: linux -*- */

/* tests/test-sparse0.C	(Formerly test-sparse-matrix.C)
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 * 2002-04-03: William J. Turner <wjturner@acm.org>
 *
 * changed name of sparse-matrix file.
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <strstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/sparse0.h"
#include "linbox/vector/random.h"
#include "linbox/util/vector-factory.h"
#include "linbox/field/vector-domain.h"
#include "linbox/randiter/nonzero.h"

#include "test-common.h"

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

template <class Field, class Vector, class Row>
static bool testIdentityApply (Field &F, const char *text, VectorFactory<Vector> &factory) 
{
	typedef SparseMatrix0 <Field, Row, Vector> Blackbox;

	char buf[80];
	strstream str (buf, 80);
	str << "Testing identity apply (" << text << ")" << ends;
	commentator.start (buf, "testIdentityApply", factory.m ());

	bool ret = true;
	bool iter_passed = true;

	VectorDomain<Field> VD (F);
	Blackbox A (F, factory.n (), factory.n ());
	Vector v, w;

	VectorWrapper::ensureDim (v, factory.n ());
	VectorWrapper::ensureDim (w, factory.n ());

	size_t i;
	typename Field::Element e;

	F.init (e, 1);

	for (i = 0; i < factory.n (); i++)
		A.put_value (pair<size_t, size_t> (i, i), e);

	while (factory) {
		commentator.startIteration (factory.j ());

		iter_passed = true;

		factory.next (v);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		printVector<Field> (F, report, v);

		A.apply (w, v);

		commentator.indent (report);
		report << "Output vector: ";
		printVector<Field> (F, report, w);

		if (!VD.areEqual (v, w))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

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

template <class Field, class Vector, class Row>
static bool testNilpotentApply (Field &F, const char *text, VectorFactory<Vector> &factory) 
{
	typedef SparseMatrix0 <Field, Row, Vector> Blackbox;

	char buf[80];
	strstream str (buf, 80);
	str << "Testing nilpotent apply (" << text << ")" << ends;
	commentator.start (buf, "testNilpotentApply", factory.m ());

	bool ret = true;
	bool iter_passed;
	Blackbox A (F, factory.n (), factory.n ());

	size_t i, j;
	typename Field::Element e;
	F.init (e, 1);
	bool even = false;

	for (i = 1; i < factory.n (); i++)
		A.put_value (pair<size_t, size_t> (i - 1, i), e);

	VectorDomain<Field> VD (F);
	Vector v, w;

	VectorWrapper::ensureDim (v, factory.n ());
	VectorWrapper::ensureDim (w, factory.n ());

	while (factory) {
		commentator.startIteration (factory.j ());

		iter_passed = true;
		even = false;

		factory.next (v);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		printVector<Field> (F, report, v);

		for (j = 0; j < factory.n () - 1; j++, even = !even)
			if (even)
				A.apply (v, w);
			else
				A.apply (w, v);

		commentator.indent (report);
		report << "A^(n-1) v:     ";
		printVector<Field> (F, report, even ? w : v);

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
		printVector<Field> (F, report, even ? v : w);

		ret = iter_passed = allZero(F, even ? v : w);

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: A^n v is non-zero" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

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

template <class Field>
bool testRandomApply1 (Field &F, size_t n, size_t iterations, size_t K) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix0 <Field, Row, Vector> Blackbox;

	commentator.start ("Testing sparse random apply (1)", "testRandomApply1", iterations);

	bool ret = true;
	bool iter_passed;

	size_t i, j, k;

	typename Field::RandIter r (F);
	typename Field::Element x;

	integer c;
	long width;

	F.characteristic (c);
	width = logp (c, 10) + 1;

	if (K > n) K = n;

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

		iter_passed = true;

		Blackbox A (F, n, n);

		for (j = 0; j < n; j++) {
			for (k = 0; k < K; k++) {
				pair<size_t, size_t> p (j, 0);

				do
					p.second = rand () % n;
				while (!F.isZero (A[p]));

				r.random (x);
				A.put_value (p, x);
			}
		}

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.prettyPrint (report, 6, width);

		for (j = 0; j < n; j++) {
			Vector e_j(n), w(n);

			F.init (e_j[j], 1);

			A.apply (w, e_j);

			for (k = 0; k < n; k++)
				if (!F.areEqual (A[pair<size_t, size_t>(k, j)], w[k]))
					ret = iter_passed = false;

			commentator.indent (report);
			report << "Output vector " << j << ": ";
			printVector<Field> (F, report, w);
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

template <class Field>
bool testRandomApply2 (Field &F, size_t n, size_t iterations, size_t N) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix0 <Field, Row, Vector> Blackbox;

	commentator.start ("Testing sparse random apply (2)", "testRandomApply2", iterations);

	bool ret = true;
	bool iter_passed;

	size_t i, j, k;

	typename Field::RandIter r (F);
	typename Field::Element x;

	integer c;
	long width;

	F.characteristic (c);
	width = logp (c, 10) + 1;

	if (N > n * n) N = n * n;

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

		iter_passed = true;

		Blackbox A (F, n, n);

		for (k = 0; k < N; k++) {
			pair<size_t, size_t> p (j, 0);

			do {
				p.first = rand () % n;
				p.second = rand () % n;
			} while (!F.isZero (A[p]));

			r.random (x);
			A.put_value (p, x);
		}

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.prettyPrint (report, 6, width);

		for (j = 0; j < n; j++) {
			Vector e_j(n), w(n);

			F.init (e_j[j], 1);

			A.apply (w, e_j);

			for (k = 0; k < n; k++)
				if (!F.areEqual (A[pair<size_t, size_t>(k, j)], w[k]))
					ret = iter_passed = false;

			commentator.indent (report);
			report << "Output vector " << j << ": ";
			printVector<Field> (F, report, w);
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

template <class Field>
bool testRandomApply3 (Field &F, size_t n, size_t iterations, size_t K) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix0 <Field, Row, Vector> Blackbox;

	commentator.start ("Testing sparse random apply (3)", "testRandomApply3", iterations);

	bool ret = true;
	bool iter_passed;

	size_t i, j, k;

	typename Field::RandIter r (F);
	typename Field::Element x, sum;

	integer c;
	long width;

	F.characteristic (c);
	width = logp (c, 10) + 1;

	if (K > n) K = n;

	Vector v(n), w(n);

	for (k = 0; k < n; k++)
		F.init (v[k], 1);

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

		iter_passed = true;

		Blackbox A (F, n, n);

		for (j = 0; j < n; j++) {
			for (k = 0; k < K; k++) {
				pair<size_t, size_t> p (j, 0);

				do
					p.second = rand () % n;
				while (!F.isZero (A[p]));

				r.random (x);
				A.put_value (p, x);
			}
		}

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.prettyPrint (report, 6, width);

		A.apply (w, v);

		for (j = 0; j < n; j++) {
			F.init (sum, 0);

			for (k = 0; k < n; k++)
				F.addin (sum, A[pair<size_t, size_t>(j, k)]);

			if (!F.areEqual (sum, w[j]))
				ret = iter_passed = false;
		}

		commentator.indent (report);
		report << "Output vector: ";
		printVector<Field> (F, report, w);

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Output vector was incorrect" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomApply3");

	return ret;
}

template <class Field> struct comp_w_ind
{ 
        bool operator() (const std::pair< size_t, typename Field::Element >& entry, 
                        size_t col_in)
        { return entry.first < col_in; }
};

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

	typedef	Modular<long>	Field;
	typedef Field::Element	Element;

	typedef std::vector<Element> Vector1;
	typedef std::vector<std::pair <size_t, Element> > Vector2;
	typedef std::map <size_t, Element > Vector3;

	typedef std::list <pair <size_t, Element> > Row;

	parseArguments (argc, argv, args);
	Modular<long> F (q);

#if 0	// The following is from the old module to test to get this one to 
	// work.  It is not intended to be kept forever.  wjt

        LinBox::SparseMatrix0Base <Field, Row> 
                B(*LinBox::newSparsemat<Field, Row>(F));

        size_t m; 
        m = B.get_rowdim();
        n = B.get_coldim();

        cout << "The sparesemat matrix contains the following Elements:" 
                << endl << B;

        cout << "Enter a vector to be multiplied by the SparseMatrix0 matrix." << endl
                << "Input the vector by entering index and value." << endl
                << "Remember matrices and vectors are indexed starting at 0." << endl
                << "End with a index of -1." << endl;

        Element zero;
        F.init(zero, 0);
        Element elem(zero);

        Vector1 x1(n, zero), y1(m, zero);
        Vector2 x2, y2;
	Vector3 x3, y3;

        Vector2::iterator x2_iter;
        Vector3::iterator x3_iter;

        bool found;

        size_t i;

        while (cin >> i)
        {
                // return also if row index is not positive
                if(i == size_t(-1)) break; 
                
                F.read(cin, elem);

                // Record element in dense vector
                x1[i] = elem;

                // Record element in sparse sequence vector

                // find appropriate location of element
                if( x2.begin() == x2.end() )
                        x2_iter = x2.end();
                else
                        x2_iter = lower_bound( x2.begin(), x2.end(), i, comp_w_ind<Field>() );

                // Check to see if element already exists.
                if ( x2.end() == x2_iter )
                        found = false;
                else if ( x2_iter->first != i )
                        found = false;
                else 
                        found = true;

                // If element is already in row, replace old value with new.
                // Otherwise, insert the element in the row.
                if (found) 
                {
                        if (F.isZero(elem))
                                x2.erase(x2_iter);
                        else
                                x2_iter->second = elem;
                } // if (found)
                else if (!F.isZero(elem))
                        x2.insert(x2_iter, make_pair(i,elem));

                // Record element in sparse associative vector

                // Find element in associative container.  
                // If exists, replace value if not zero, or remove if value is zero.
                // If not found, insert non-zero element
                if ( (x3_iter = x3.find(i)) != x3.end() )
                {
                        if (F.isZero(elem))
                                x3.erase(x3_iter);
                        else
                                x3_iter->second = elem;
                }
                else
                {
                        if (!F.isZero(elem))
                        x3.insert(make_pair(i, elem));
                }

        } // while (cin >> i)

        cout << "*** Running tests with dense vector." << endl;

        cout << "Dense vector x1:" << endl;
	printVector(F, cout, x1);


        LinBox::SparseMatrix0<Field, Row, Vector1> S1(B);
        LinBox::BlackboxArchetype<Vector1>& A1 = S1;

        A1.apply(y1,x1);

        cout << "Using A1.apply(y1,x1) gives vector y1:" << endl;
	printVector(F, cout, y1);

        y1 = Vector1(m, zero);
        y1 = A1.apply(x1);

        cout << "Using y1 = A1.apply(x1) gives vector y1:" << endl;
	printVector(F, cout, y1);

        y1 = x1;
        y1 = A1.applyIn(y1);

        cout << "Using y1 = x1; y1 = A1.applyIn(y1) gives vector y1:" << endl;
	printVector(F, cout, y1);

        y1 = Vector1(m, zero);
        A1.applyTranspose(y1,x1);

        cout << "Using A1.applyTranspose(y1,x1) gives vector y1:" << endl;
	printVector(F, cout, y1);

        y1 = Vector1(m, zero);
        y1 = A1.applyTranspose(x1);

        cout << "Using y1 = A1.applyTranspose(x1) gives vector y1:" << endl;
	printVector(F, cout, y1);

        y1 = x1;
        y1 = A1.applyTransposeIn(y1);

        cout << "Using y1 = x1; y1 = A1.applyTransposeIn(y1) gives vector y1:" << endl;
	printVector(F, cout, y1);

        cout << "*** Running tests with sparse sequence vector." << endl;

        cout << "Sparse sequence vector x2:" << endl;
	printVector(F, cout, x2);

        LinBox::SparseMatrix0<Field, Row, Vector2> S2(B);
        LinBox::BlackboxArchetype<Vector2>& A2 = S2;
        
        A2.apply(y2,x2);

        cout << "Using A2.apply(y2,x2) gives vector y2:" << endl;
	printVector(F, cout, y2);

        y2 = Vector2();
        y2 = A2.apply(x2);

        cout << "Using y2 = A2.apply(x2) gives vector y2:" << endl;
	printVector(F, cout, y2);

        y2 = x2;
        y2 = A2.applyIn(y2);

        cout << "Using y2 = x2; y2 = A2.applyIn(y2) gives vector y2:" << endl;
	printVector(F, cout, y2);

        y2 = Vector2();
        A2.applyTranspose(y2,x2);

        cout << "Using A2.applyTranspose(y2,x2) gives vector y2:" << endl;
	printVector(F, cout, y2);

        y2 = Vector2();
        y2 = A2.applyTranspose(x2);

        cout << "Using y2 = A2.applyTranspose(x2) gives vector y2:" << endl;
	printVector(F, cout, y2);

        y2 = x2;
        y2 = A2.applyTransposeIn(y2);

        cout << "Using y2 = x2; y2 = A2.applyTransposeIn(y2) gives vector y2:" << endl;
	printVector(F, cout, y2);

	cout << "*** Running tests with sparse associative vector." << endl;

	cout << "Sparse associative vector x3:" << endl;
	printVector(F, cout, x3);

	LinBox::SparseMatrix0<Field, Row, Vector3> S3(B);
	LinBox::BlackboxArchetype<Vector3>& A3 = S3;
	
	A3.apply(y3,x3);

	cout << "Using A3.apply(y3,x3) gives vector y3:" << endl;
	printVector(F, cout, y3);

	y3 = Vector3();
	y3 = A3.apply(x3);

	cout << "Using y3 = A3.apply(x3) gives vector y3:" << endl;
	printVector(F, cout, y3);

	y3 = x3;
	y3 = A3.applyIn(y3);

	cout << "Using y3 = x3; y3 = A3.applyIn(y3) gives vector y3:" << endl;
	printVector(F, cout, y3);

	y3 = Vector3();
	A3.applyTranspose(y3,x3);

	cout << "Using A3.applyTranspose(y3,x3) gives vector y3:" << endl;
	printVector(F, cout, y3);

	y3 = Vector3();
	y3 = A3.applyTranspose(x3);

	cout << "Using y3 = A3.applyTranspose(x3) gives vector y3:" << endl;
	printVector(F, cout, y3);

	y3 = x3;
	y3 = A3.applyTransposeIn(y3);

	cout << "Using y3 = x3; y3 = A3.applyTransposeIn(y3) gives vector y3:" << endl;
	printVector(F, cout, y3);

#endif // end of stuff from old module

	srand (time (NULL));

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	cout << "Sparse matrix black box test suite" << endl << endl;

	RandomDenseVectorFactory<Field> factory1 (F, n, iterations);
	RandomSparseSeqVectorFactory<Field> factory2 (F, n, n / 10, iterations);
	RandomSparseMapVectorFactory<Field> factory3 (F, n, n / 10, iterations);

	RandomDenseVectorFactory<Field, NonzeroRandIter<Field> > factory4 (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), n, iterations);
	RandomSparseSeqVectorFactory<Field, NonzeroRandIter<Field> > factory5 (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), n, n / 10, iterations);
	RandomSparseMapVectorFactory<Field, NonzeroRandIter<Field> > factory6 (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), n, n / 10, iterations);

	cout << "Running for dense vectors" << endl;
	if (!testIdentityApply<Field, Vector1, Row>  (F, "sparse sequence/dense", factory1)) pass = false;
	if (!testNilpotentApply<Field, Vector1, Row> (F, "sparse sequence/dense", factory4)) pass = false;
	if (!testRandomApply1<Modular<long> >     (F, n, iterations, k)) pass = false;
	if (!testRandomApply2<Modular<long> >     (F, n, iterations, N)) pass = false;
	if (!testRandomApply3<Modular<long> >     (F, n, iterations, k)) pass = false;

	cout << "Running for sparse sequence vectors" << endl;
	if (!testIdentityApply<Field, Vector2, Row>  (F, "sparse sequence/sparse sequence", factory2)) pass = false;
	if (!testNilpotentApply<Field, Vector2, Row> (F, "sparse sequence/sparse sequence", factory5)) pass = false;
	if (!testRandomApply1<Modular<long> >     (F, n, iterations, k)) pass = false;
	if (!testRandomApply2<Modular<long> >     (F, n, iterations, N)) pass = false;
	if (!testRandomApply3<Modular<long> >     (F, n, iterations, k)) pass = false;

//	if (!testIdentityApply<Field, Vector3, Row>    (F, n, iterations)) pass = false;
//	if (!testIdentityApply<Field, Vector4, Row>    (F, n, iterations)) pass = false;
	

	return pass ? 0 : -1;
}
