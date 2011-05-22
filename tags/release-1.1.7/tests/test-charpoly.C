/* tests/test-charpoly.C
 * Copyright (C) LinBox
 * Written by bds (starting from test-charpoly.C)
 *
 * See COPYING for license information.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdio>
template<class T, template <class T> class Container>
std::ostream& operator<< (std::ostream& o, const Container<T>& C) {
	for(typename Container<T>::const_iterator refs =  C.begin();
	    refs != C.end() ;
	    ++refs )
		o << (*refs) << " " ;
	return o << std::endl;
}

#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/solutions/charpoly.h"
#include "linbox/util/commentator.h"
//#include "linbox/ring/givaro-polynomial.h"
#include "linbox/vector/stream.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: charpoly of the identity matrix
 *
 * Construct the identity matrix and compute its characteristic polynomial. Confirm
 * that the result is (x-1)^n
 *
 * Z - Integers representation over which to perform computations
 * n - Dimension to which to make matrix
 *
 * Return true on success and false on failure
 */

template <class Dom, class Polynomial>
typename Dom::Element eval (const Dom& D,
			    typename Dom::Element& value,
			    const Polynomial& P,
			    const typename Dom::Element x){
	typename Dom::Element tmp = P[P.size()-1];
	for (int i = P.size()-2; i >= 0; --i){
		D.mulin (tmp, x);
		D.addin (tmp, P[i]);
	}
	return value = tmp; 
}

template <class Dom>
static bool testIdentityCharpoly (Dom &Z, size_t n, bool symmetrizing=false) 
{
	typedef typename Dom::Element Element;
	typedef vector<Element> Vector;
	typedef ScalarMatrix<Dom> Blackbox;
	typedef GivPolynomialRing<Dom, Dense> PolDom;
	typedef typename PolDom::Element Polynomial;
	//typedef Vector Polynomial;

	commentator.start ("Testing identity Charpoly", "testIdentityCharpoly");

	bool ret = true;
	Element one; Z.init(one, 1);
	Element negone; Z.init(negone, -1);

	//PolDom IPD(Z);

	Blackbox A (Z, n, one);

	Polynomial phi;

	charpoly (phi, A);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Characteristic polynomial is: ";
	printPolynomial<Dom, Polynomial> (Z, report, phi);

	// partial check - just that charpoly has right values at 0, 1, -1.
	Element val, val2, neg2, pow2;
	// value at 1 should be zero
	eval(Z, val, phi, one);
	if (! Z.isZero(val) ) ret = false;
	// value at zero should be (-1)^n
	val = (n % 2 == 0) ? one : negone;
	if (! Z.areEqual(val, phi[0])) ret = false;
	// value at -1 should be (-2)^n
	eval(Z, val2, phi, negone);
	Z.init(neg2, -2); Z.init(pow2, 1);
	for (size_t i = 0; i < n; ++i) Z.mulin(pow2, neg2);
	if (! Z.areEqual(val2, pow2)) ret = false;

	if (! ret){
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Characteristic polynomial is incorrect" << endl;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testIdentityCharpoly");

	return ret;
}

/* Test 2: characteristic polynomial of a nilpotent matrix
 *
 * Construct an index-n nilpotent matrix and compute its characteristic
 * polynomial. Confirm that the result is x^n
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testNilpotentCharpoly (Field &F, size_t n)
{
	typedef vector <typename Field::Element> Vector;
// 	typedef GivPolynomialRing<Field, Dense> PolDom;
// 	typedef typename PolDom::Element Polynomial;
	typedef Vector Polynomial;
	typedef pair <vector <size_t>, vector <typename Field::Element> > Row;
	typedef SparseMatrix <Field> Blackbox;

	commentator.start ("Testing nilpotent charpoly", "testNilpotentCharpoly");

	bool ret = true;
	bool lowerTermsCorrect = true;

	size_t i;

	StandardBasisStream<Field, Row> stream (F, n);
	Row v;
	stream.next (v);
	Blackbox A (F, stream);

	Polynomial phi;

	charpoly (phi, A);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "characteristic polynomial is: ";
	printPolynomial (F, report, phi);

	for (i = 0; i < n - 1; i++)
		if (!F.isZero (phi[i]))
			lowerTermsCorrect = false;

	if (phi.size () != n + 1 || !F.isOne (phi[n]) || !lowerTermsCorrect) {
		ret = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: characteristic polynomial is incorrect" << endl;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testNilpotentCharpoly");

	return ret;
}

#if 1
/* Test 3: Random charpoly of sparse matrix
 *
 * Generates a random sparse matrix with K nonzero elements per row and computes
 * its characteristic polynomial. Then computes random vectors and applies the
 * polynomial to them in Horner style, checking whether the result is 0.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * K - Number of nonzero elements per row
 * numVectors - Number of random vectors to which to apply the characteristic polynomial
 *
 * Return true on success and false on failure
 */

template <class Field, class Row, class Vector>
bool testRandomCharpoly (Field                 &F,
			int                    iterations,
			VectorStream<Row>    &A_stream,
			VectorStream<Vector> &v_stream)
{
	//typedef GivPolynomialRing<Field, Dense> PolDom;
	//typedef typename PolDom::Element Polynomial;
	typedef std::vector<typename Field::Element> Polynomial;
	typedef SparseMatrix <Field> Blackbox;

	commentator.start ("Testing sparse random charpoly", "testRandomCharpoly", iterations);

	bool ret = true;
	bool iter_passed;

	VectorDomain<Field> VD (F);

	Vector v, w;

	VectorWrapper::ensureDim (v, v_stream.n ());
	VectorWrapper::ensureDim (w, v_stream.n ());

	for (int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		iter_passed = true;

		A_stream.reset ();
		Blackbox A (F, A_stream);

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.write (report, FORMAT_PRETTY);

		Polynomial phi;

		charpoly (phi, A);

		report << "characteristic polynomial is: ";
		printPolynomial (F, report, phi);

		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "deg charpoly (A) = " << phi.size () - 1 << endl;

		v_stream.reset ();

		while (v_stream) {
			v_stream.next (v);

			report << "Input vector  " << v_stream.j () << ": ";
			VD.write (report, v);
			report << endl;

			applyPoly (F, w, A, phi, v);

			report << "Output vector " << v_stream.j () << ": ";
			VD.write (report, w);
			report << endl;

			if (!VD.isZero (w))
				ret = iter_passed = false;
		}

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Output vector was incorrect" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomCharpoly");

	return ret;
}
#endif

int main (int argc, char **argv)
{
	bool pass = true;

	std::cout<<setprecision(8);
	std::cerr<<setprecision(8);
	static size_t n = 50;
	static integer q = 33554467U; 
	//static integer q = 1000003U; // 33554467U; 
	static int iterations = 1;
	static int numVectors = 100;
	static int k = 3;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.",                 TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].",          TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.",                    TYPE_INT,     &iterations },
		{ 'v', "-v V", "Use V test vectors for the random charpoly tests.",      TYPE_INT,     &numVectors },
		{ 'k', "-k K", "K nonzero Elements per row in sparse random apply test.", TYPE_INT,     &k },
		{ '\0' }
	};


	parseArguments (argc, argv, args);

	// Temporarily, only Modular<double> is enabled for the givaro/ntl factorization based charpoly
	typedef Modular<double> Field;
	typedef vector<Field::Element> DenseVector;
	typedef SparseMatrix<Field>::Row SparseVector;
	//typedef pair<vector<size_t>, vector<Field::Element> > SparseVector;
	Field F (q);
	srand (time (NULL));

	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << endl << "Black box characteristic polynomial of a matrix over a prime field test suite" << endl;

	RandomDenseStream<Field, DenseVector, NonzeroRandIter<Field> >
		v_stream (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), n, numVectors);
	RandomSparseStream<Field, SparseVector, NonzeroRandIter<Field> >
		A_stream (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), (double) k / (double) n, n, n);

	if (!testNilpotentCharpoly (F, n)) pass = false;
	if (!testRandomCharpoly    (F, iterations, A_stream, v_stream)) pass = false;

	// symmetrizing
	if (!testIdentityCharpoly  (F, n, true)) pass = false;
	//need other tests...

	typedef vector<PID_integer::Element> ZDenseVector;
	typedef SparseMatrix<PID_integer>::Row ZSparseVector;
	typedef pair<vector<size_t>, vector<Field::Element> > SparseVector;
	PID_integer  Z;
	srand (time (NULL));

	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	report << endl << "Black box characteristic polynomial of an integer matrix test suite" << endl;

	RandomDenseStream<PID_integer, ZDenseVector, NonzeroRandIter<PID_integer> >
		zv_stream (Z, NonzeroRandIter<PID_integer> (Z, PID_integer::RandIter (Z)), n, numVectors);
	RandomSparseStream<PID_integer, SparseVector, NonzeroRandIter<PID_integer> >
		zA_stream (Z, NonzeroRandIter<PID_integer> (Z, PID_integer::RandIter (Z)), (double) k / (double) n, n, n);

	//no symmetrizing
	if (!testIdentityCharpoly  (Z, n)) pass = false;
	if (!testNilpotentCharpoly (Z, n)) pass = false;

	//Comment by Z. Wan. Stream doesn't work here
	//if (!testRandomCharpoly    (Z, iterations, zA_stream, zv_stream)) pass = false;

	// symmetrizing
	if (!testIdentityCharpoly  (Z, n, true)) pass = false;
	//need other tests...

	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
