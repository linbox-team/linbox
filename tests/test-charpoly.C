/* tests/test-charpoly.C
 * Copyright (C) LinBox
 * Written by bds ++
 *
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file  tests/test-charpoly.C
 * @ingroup tests
 * @brief tests the characteristic polynomial of sparse and special matrices
 * @warning gcc-4.2 produces bad optimized code there
 * @bug occasionnnaly there is a "SIGFPE, Arithmetic exception." in CRA
 * @bug testRandomCharpoly is not always tested !!
 * @test characteristic polynomial of some matrices (sparse, special)
 */


#include "linbox/linbox-config.h"

#include <iostream>
#include <iomanip>
#include <fstream>

#include <cstdio>

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/solutions/charpoly.h"
#include "linbox/util/commentator.h"
#include "linbox/ring/polynomial-ring.h"
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
			    const typename Dom::Element x)
{
	typename Dom::Element tmp = P[P.size()-1];
	for (int i = (int)P.size()-2; i >= 0; --i){
		D.mulin (tmp, x);
		D.addin (tmp, P[(size_t)i]);
	}
	return value = tmp;
}

template <class Dom>
static bool testIdentityCharpoly (Dom &Z, size_t n, bool symmetrizing=false)
{
	typedef typename Dom::Element Element;
	typedef ScalarMatrix<Dom> Blackbox;
//	typedef GivPolynomialRing<Dom, Givaro::Dense> PolDom;
//  typedef BlasVector<Dom,GivPolynomialRing<Dom, Givaro::Dense> > PolDom;
    typedef DensePolynomial<Dom> Polynomial;

	LinBox::commentator().start ("Testing identity Charpoly", "testIdentityCharpoly");

	bool ret = true;

	//PolDom IPD(Z);

	Blackbox A (Z, n, n, Z.one);

	Polynomial phi(Z);

	charpoly (phi, A);

	ostream &report = LinBox::commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Characteristic polynomial is: ";
	printPolynomial<Dom, Polynomial> (Z, report, phi);

	// partial check - just that charpoly has right values at 0, 1, -1.
	Element val, val2, neg2, pow2;
	// value at 1 should be zero
	eval(Z, val, phi, Z.one);
	if (! Z.isZero(val) ) ret = false;
	// value at zero should be (-1)^n
	val = (n % 2 == 0) ? Z.one : Z.mOne;
	if (! Z.areEqual(val, phi[0])) ret = false;
	// value at -1 should be (-2)^n
	eval(Z, val2, phi, Z.mOne);
	Z.init(neg2, -2_i64); Z.assign(pow2, Z.one);
	for (size_t i = 0; i < n; ++i) Z.mulin(pow2, neg2);
	if (! Z.areEqual(val2, pow2)) ret = false;

	if (! ret){
		LinBox::commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Characteristic polynomial is incorrect" << endl;
	}

	LinBox::commentator().stop (MSG_STATUS (ret), (const char *) 0, "testIdentityCharpoly");

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
//	typedef GivPolynomialRing<Field, Givaro::Dense> PolDom;
//	typedef typename PolDom::Element Polynomial;
	typedef DensePolynomial<Field> Polynomial;
	typedef std::pair <std::vector <size_t>, std::vector <typename Field::Element> > Row;
	typedef SparseMatrix<Field, typename VectorTraits<Row>::SparseFormat> Blackbox;

	LinBox::commentator().start ("Testing nilpotent charpoly", "testNilpotentCharpoly");

	bool ret = true;
	bool lowerTermsCorrect = true;

	size_t i;

	StandardBasisStream<Field, Row> stream (F, n);
	Row v;
	stream.next (v);
	Blackbox A (F, stream);

	ostream &who = LinBox::commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	who << "Matrix:" << endl;
	A.write (who, Tag::FileFormat::Pretty);

        Polynomial phi(F);

	charpoly (phi, A);

	ostream &report = LinBox::commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "characteristic polynomial is: ";
	printPolynomial (F, report, phi);

	linbox_check(n);
	for (i = 0; i < n - 1; i++)
		if (!F.isZero (phi[(size_t)i]))
			lowerTermsCorrect = false;

	if (phi.size () != n + 1 || !F.isOne (phi[n]) || !lowerTermsCorrect) {
		ret = false;
		LinBox::commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: characteristic polynomial is incorrect (should be x^" << n << ')' << endl;
	}

	LinBox::commentator().stop (MSG_STATUS (ret), (const char *) 0, "testNilpotentCharpoly");

	return ret;
}

static bool testSageBug(){

        Givaro::ZRing<Givaro::Integer> Z;
        DenseMatrix<Givaro::ZRing<Givaro::Integer> > A(Z,4,4);
        for (uint32_t i=0; i<4; ++i)
                for (uint32_t j=0; j<4; ++j)
                        A.setEntry(i,j, Givaro::Integer(i*4+j+1));
        typedef DensePolynomial<Givaro::ZRing<Givaro::Integer> > Polynomial;
        Polynomial phi(Z);
        charpoly(phi,A);
        if (Z.areEqual(phi[0],0) &&
            Z.areEqual(phi[1],0) &&
            Z.areEqual(phi[2],-80) &&
            Z.areEqual(phi[3],-34) &&
            Z.areEqual(phi[4],1) )
            return true;
        else return false;
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
			VectorStream<Row>    &A_stream,
			VectorStream<Vector> &v_stream)
{
//     typedef GivPolynomialRing<Field, Givaro::Dense> PolDom;
        typedef DensePolynomial<Field> Polynomial;
	typedef SparseMatrix<Field> Blackbox;

	LinBox::commentator().start ("Testing sparse random charpoly", "testRandomCharpoly", 1);

	bool ret = true;

	VectorDomain<Field> VD (F);
	Vector w(F), v(F);
	VectorWrapper::ensureDim (v, v_stream.n ());
	VectorWrapper::ensureDim (w, v_stream.n ());

	A_stream.reset ();
	Blackbox A (F, A_stream);

	ostream &report = LinBox::commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
	report << "Matrix:" << endl;
	A.write (report, Tag::FileFormat::Pretty);

    Polynomial phi(F);

	charpoly (phi, A);

	report << "characteristic polynomial is: ";
	printPolynomial (F, report, phi);

	LinBox::commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
		<< "deg charpoly (A) = " << phi.size () - 1 << endl;

	v_stream.reset ();

	while (v_stream) {
		v_stream.next (v);
		VD.write (report << "Input vector  " << v_stream.j () << ": ", v) << endl;

		applyPoly (F, w, A, phi, v);
		VD.write (report << "Output vector " << v_stream.j () << ": ", w) << endl;

		if (!VD.isZero (w)) { ret = false; break; }
	}

	if (!ret)
		LinBox::commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Output vector was incorrect" << endl;

    typedef PolynomialRing<Field> PolyDom;
    Polynomial psi(F);

	charpoly (psi, A, Method::Blackbox() );

    ret = ret && PolyDom(F).areEqual(phi, psi);

	if (!ret) {
		LinBox::commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Auto charpoly differs from BB charpoly" << endl;

        PolyDom(F).write(report << "phi: ", phi) << std::endl;
        PolyDom(F).write(report << "psi: ", psi) << std::endl;
    }


	LinBox::commentator().stop (MSG_STATUS (ret), (const char *) 0, "testRandomCharpoly");
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
	//static integer q = 103U; // 33554467U;
	static int numVectors = 10;
	static int k = 3;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.",                 TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].",          TYPE_INTEGER, &q },
		{ 'v', "-v V", "Use V test vectors for the random charpoly tests.",      TYPE_INT,     &numVectors },
		{ 'k', "-k K", "K nonzero Elements per row in sparse random apply test.", TYPE_INT,     &k },
		END_OF_ARGUMENTS
	};


	parseArguments (argc, argv, args);

	LinBox::commentator().getMessageClass (TIMING_MEASURE).setMaxDepth (10);
	LinBox::commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
	LinBox::commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	ostream &report = LinBox::commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	srand ((unsigned)time (NULL));

/**************/
	report << endl << "Black box characteristic polynomial over a prime field test suite" << endl;

	// Temporarily, only Givaro::Modular<double> is enabled for the givaro/ntl factorization based charpoly
	typedef Givaro::Modular<double> Field;
	typedef BlasVector<Field> DenseVector;
	typedef SparseMatrix<Field>::Row SparseVector;
	Field F (q);
	Field::RandIter mygen(F);
        Givaro::GeneralRingNonZeroRandIter<Field> myNZgen(mygen);

	RandomDenseStream<Field, DenseVector, Givaro::GeneralRingNonZeroRandIter<Field> >
        v_stream (F, myNZgen, n, (size_t)numVectors);
	RandomSparseStream<Field, SparseVector, Givaro::GeneralRingNonZeroRandIter<Field> >
	A_stream (F, myNZgen, (double) k / (double) n, n, n);

	if (!testNilpotentCharpoly (F, n)) pass = false;
	if (!testRandomCharpoly    (F, A_stream, v_stream)) pass = false;

	// symmetrizing
	if (!testIdentityCharpoly  (F, n, true)) pass = false;
	//need other tests...

/**************/
	report << endl << "Black box characteristic polynomial over Z test suite" << endl;

	Givaro::ZRing<integer>  Z ;
	typedef BlasVector<Givaro::ZRing<integer>> ZDenseVector;

	Givaro::ZRing<integer>::RandIter myZgen(Z);
        Givaro::GeneralRingNonZeroRandIter<Givaro::ZRing<integer>> myNzZgen(myZgen);

	RandomDenseStream<Givaro::ZRing<integer>, ZDenseVector, Givaro::GeneralRingNonZeroRandIter<Givaro::ZRing<integer>> >
        zv_stream (Z, myNzZgen, n, (size_t)numVectors);
	RandomSparseStream<Givaro::ZRing<integer>, SparseVector, Givaro::GeneralRingNonZeroRandIter<Givaro::ZRing<integer>> >
	zA_stream (Z, myNzZgen, (double) k / (double) n, n, n);

	//no symmetrizing
	if (!testIdentityCharpoly  (Z, n)) pass = false;
	if (!testNilpotentCharpoly (Z, n)) pass = false;

	//Comment by Z. Wan. Stream doesn't work here
//	if (!testRandomCharpoly    (Z, zA_stream, zv_stream)) pass = false;

	// symmetrizing
	if (!testIdentityCharpoly  (Z, n, true)) pass = false;
	//need other tests...

	if (not testSageBug()) pass = false;

	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
