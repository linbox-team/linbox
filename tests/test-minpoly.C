/* Copyright (C) the LinBox group
 * tests/test-minpoly.C *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 * 2002-04-03: William J. Turner <wjturner@acm.org>
 *
 * changed name of sparse-matrix file.
 *
 * --------------------------------------------------------
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


/*! @file  tests/test-minpoly.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */



#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>

#include <cstdio>

#include "givaro/modular.h"
#include "givaro/gfq.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/polynomial/dense-polynomial.h"
#include "linbox/util/commentator.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/solutions/charpoly.h"
#include "linbox/vector/stream.h"

#include "linbox/vector/blas-vector.h"
#include "fflas-ffpack/utils/test-utils.h"
#include "test-common.h"

using namespace LinBox;

/* Test 0: Minimal polynomial of the zero matrix
 */
template <class Field, class Meth>
static bool testZeroMinpoly (Field &F, size_t n, const Meth& M)
{
	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	commentator().start ("Testing zero minpoly", "testZeroMinpoly");
//	typedef BlasVector<Field> Polynomial;
	typedef DensePolynomial<Field> Polynomial;
	Polynomial phi(F);
	SparseMatrix<Field> A(F, n, n);
	minpoly(phi, A, M);

    A.write(report, Tag::FileFormat::Maple);
	report << "Minimal polynomial is: ";

	printPolynomial<Field, Polynomial> (F, report, phi);

	bool ret;
	if (phi.size () == 2 && F.isZero (phi[0]) && F.isOne(phi[1]) )
		ret = true;
	else
	{
		report << "ERROR: A = 0, should get x, got ";
		printPolynomial<Field, Polynomial> (F, report, phi);
		report << endl;
		ret = false;
	}
	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testZeroMinpoly");
	return ret;
}

/* Test 1: Minimal polynomial of the identity matrix
 *
 * Construct the identity matrix and compute its minimal polynomial. Confirm
 * that the result is x-1
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 *
 * Return true on success and false on failure
 */

template <class Field, class Meth>
static bool testIdentityMinpoly (Field &F, size_t n, const Meth& M)
{
//	typedef BlasVector<Field> Polynomial;
	typedef DensePolynomial<Field> Polynomial;
	typedef ScalarMatrix<Field> Blackbox;

	commentator().start ("Testing identity minpoly", "testIdentityMinpoly");

	typename Field::Element c0, c1;

        //StandardBasisStream<Field, Row> stream (F, n);
	Blackbox A (F, n, n, F.one);

	Polynomial phi(F);

	minpoly (phi, A, M );

	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Minimal polynomial is: ";
	printPolynomial<Field, Polynomial> (F, report, phi);

	bool ret;

	F.assign(c0, F.mOne);
	F.assign(c1, F.one);

	if (phi.size () == 2 && F.areEqual (phi[0], c0) && F.areEqual (phi[1], c1))
		ret = true;
	else {
		ret = false;
		commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: A = I, should get x-1, got ";
		printPolynomial<Field, Polynomial> (F, report, phi);
	}

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testIdentityMinpoly");

	return ret;
}

/* Test 2: Minimal polynomial of a nilpotent matrix
 *
 * Construct an index-n nilpotent matrix and compute its minimal
 * polynomial. Confirm that the result is x^n
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 *
 * Return true on success and false on failure
 */

template <class Field, class Meth>
static bool testNilpotentMinpoly (Field &F, size_t n, const Meth& M)
{
//	typedef BlasVector<Field> Polynomial;
	typedef DensePolynomial<Field> Polynomial;
	typedef SparseMatrix<Field> Blackbox;
	typedef typename Blackbox::Row Row;

	commentator().start ("Testing nilpotent minpoly", "testNilpotentMinpoly");
	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	bool ret = false;
	bool lowerTermsCorrect = true;

	size_t i;

	StandardBasisStream<Field, Row> stream (F, n);
	Row v;
	stream.next (v);
	Blackbox A (F, stream); // first subdiagonal is 1's.

    report << "Matrix:" << endl;
    A.write (report, Tag::FileFormat::Maple);

	Polynomial phi(F);

	minpoly (phi, A, M);


    if (phi.size()) { // Minimal polynomial cannot be zero
        size_t d(phi.size()-1);
        for (i = 0; i < d; i++)
            if (!F.isZero (phi[i]))
                lowerTermsCorrect = false;
        if (F.isOne (phi[d]) && lowerTermsCorrect) {
                // Returned polynomial is a divisor of X^n
            ret = true;
            if (d != n) {
                    // Blackbox probabilistic methods could return only a divisor
                commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_WARNING)
                    << "Warning: A^n = 0, should get x^" << n
                    <<", got only a divisor: ";
                printPolynomial (F, report, phi);
            }
        }
    }

    if (! ret) {
        commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: A^n = 0, should get x^" << n <<", got ";
		printPolynomial (F, report, phi);
    }

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testNilpotentMinpoly");

	return ret;
}

/* Test 3: Random minpoly of sparse matrix
 *
 * Generates a random sparse matrix with K nonzero elements per row and computes
 * its minimal polynomial. Then computes random vectors and applies the
 * polynomial to them in Horner style, checking whether the result is 0.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * K - Number of nonzero elements per row
 * numVectors - Number of random vectors to which to apply the minimal polynomial
 *
 * Return true on success and false on failure
 */

template <class Field, class BBStream, class VectStream, class Meth>
bool testRandomMinpoly (Field                 &F,
                        int                    iterations,
                        BBStream    &A_stream,
                        VectStream &v_stream,
                        const Meth& M)
{
	typedef DensePolynomial<Field> Polynomial;
	typedef SparseMatrix<Field> Blackbox;
    typedef typename VectStream::Vector Vector;

    commentator().start ("Testing sparse random minpoly", "testRandomMinpoly", (unsigned int)iterations);

	bool ret = true;

	VectorDomain<Field> VD (F);

	Vector v(F), w(F);

	VectorWrapper::ensureDim (v, v_stream.n ());
	VectorWrapper::ensureDim (w, v_stream.n ());

	for (int i = 0; i < iterations; i++) {
		commentator().startIteration ((unsigned int)i);

		bool iter_passed = true;
		A_stream.reset ();
		Blackbox A (F, A_stream);

		ostream &report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.write (report, Tag::FileFormat::Maple);

		Polynomial phi(F);

		minpoly (phi, A, M );

		report << "Minimal polynomial is: ";
		printPolynomial (F, report, phi);

		commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "deg minpoly (A) = " << phi.size () - 1 << endl;

		v_stream.reset ();

		while (iter_passed && v_stream) {
			v_stream.next (v);

			report << "Input vector  " << v_stream.j () << ": ";
			VD.write (report, v);
			report << endl;

			applyPoly (F, w, A, phi, v);

			report << "Output vector " << v_stream.j () << ": ";
			VD.write (report, w);
			report << endl;

			if (!VD.isZero (w)) {
                    // Blackbox probabilistic methods could return only a divisor
                report << "polynomial is only a divisor, trying divisibility ... " << std::endl;
                Polynomial psi(F);
                Method::DenseElimination dem; dem.certifyInconsistency=true;
                charpoly (psi, A, dem );

                report << "Charpoly is: ";
                printPolynomial (F, report, psi);

                Polynomial rem(F);
                typename Polynomial::Domain_t PolDom(F);

                PolDom.mod(rem,psi,phi);
                Givaro::Degree deg; PolDom.degree(deg, rem);


                if (deg > 0) {
                    commentator().report (Commentator::LEVEL_IMPORTANT,
                                          INTERNAL_WARNING) <<
                        "ERROR: A = rand, minpoly does not divide charpoly." << endl;
                    ret = iter_passed = false;
                }
            }
		}

		if (!iter_passed)
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: A = rand, purported-minpoly(A) is not annihilator." << endl;
		commentator().stop ("done");
		commentator().progress ();
	}

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testRandomMinpoly");

	return ret;
}

/* Test 4: test gram matrix.
     A test of behaviour with self-orthogonal rows and cols.

	 Gram matrix is 1's offdiagonal and 0's on diagonal.  Self orthogonality when characteristic | n-1.

	 Arg m is ignored.
	 if p := characteristic is small then size is p+1 and minpoly is x^2 + x
	 (because A^2 = (p-1)A).
	 if p is large, size is 2 and minpoly is x^2 -1.
*/
template <class Field, class Meth>
static bool testGramMinpoly (Field &F, size_t m, const Meth& M)
{
	commentator().start ("Testing gram minpoly", "testGramMinpoly");
	typedef DensePolynomial<Field> Polynomial;
	integer p;
    size_t n;
	F.characteristic(p);
	if (p > 29u) n = 2;
	Polynomial phi(F);
	BlasMatrix<Field> A(F, n, n);
	for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) A.setEntry(i, j, F.one);
	for (size_t i = 0; i < n; ++i) A.setEntry(i, i, F.zero);
	minpoly(phi, A, M);

	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
    A.write (report);
	report << "Minimal polynomial is: ";
	printPolynomial<Field, Polynomial> (F, report, phi);

	bool ret;
	if (n == 2) {
		if ( phi.size() == 3 && F.areEqual(phi[0], F.mOne) && F.isZero(phi[1]) && F.isOne(phi[2]))
			ret = true;
		else
		{
            if ( phi.size() == 2 && F.isOne(phi[1]) &&
                 ( F.areEqual(phi[0], F.mOne) || F.isOne(phi[0]) )
                 )
            {
                commentator().report (Commentator::LEVEL_IMPORTANT,
                                      INTERNAL_WARNING) <<
                    "Warning: got only a divisor. " << endl;
                ret = true;
            } else {
                commentator().report (Commentator::LEVEL_IMPORTANT,
                                      INTERNAL_ERROR) <<
                    "ERROR: A = gram, should get x^2 - 1, got ";
                printPolynomial<Field, Polynomial> (F, report, phi);
                ret = false;
            }
		}
	} else {
		if (phi.size() == 3 && F.isZero(phi[0]) && F.isOne(phi[1]) && F.isOne(phi[2]))
			ret = true;
		else
		{
            if ( phi.size() == 2 && F.isOne(phi[1]) &&
                 ( F.isOne(phi[0]) || F.isZero(phi[0]) )
                 )
            {
                commentator().report (Commentator::LEVEL_IMPORTANT,
                                      INTERNAL_WARNING) <<
                    "Warning: got only a divisor. " << endl;
                ret = true;
            } else {
                commentator().report (Commentator::LEVEL_IMPORTANT,
                                      INTERNAL_ERROR) <<
                    "ERROR: A = gram, should get x^2 + x, got ";
                printPolynomial<Field, Polynomial> (F, report, phi);
                ret = false;
            }
		}
	}
	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testGramMinpoly");
	return ret;
}

template <class Field>
bool run_with_field(integer q, int e, uint64_t b, size_t n, int iter, int numVectors, int k, uint64_t seed){
	bool ok = true;
	int nbiter = iter;

	while (ok && nbiter)
	{
		Field* F;
        integer card=q;
        do{
            F = FFPACK::chooseField<Field>(q, b, ++seed); // F, characteristic q of b bits
            card = F->cardinality();
        } while (card < 2*uint64_t(n)*uint64_t(n) && card != 0); // ensures high probability of succes of the probabilistic algorithm
		//typename Field::RandIter G(*F, b, seed); //random generator over F
        Givaro::Integer samplesize(1); samplesize<<=b;
        typename Field::RandIter G(*F,seed,samplesize); //random generator over F
		typename Field::NonZeroRandIter NzG(G); //non-zero random generator over F

		if(F == nullptr)
			return true; //if F is null, nothing to test, just pass

        ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
        F->write(report << "Working over: ") << std::endl;
            /*
              ostringstream oss;
              F->write(oss);
              cout.fill('.');
              cout<<"Checking ";
              cout.width(40);
              cout<<oss.str();
              cout<<" ... ";
            */

		ok &= testZeroMinpoly	   (*F, n, Method::Auto());
		ok &= testZeroMinpoly	   (*F, n, Method::Elimination());
		ok &= testZeroMinpoly	   (*F, n, Method::Blackbox());
        ok &= testIdentityMinpoly  (*F, n, Method::Auto());
        ok &= testIdentityMinpoly  (*F, n, Method::Elimination());
        ok &= testIdentityMinpoly  (*F, n, Method::Blackbox());
        ok &= testNilpotentMinpoly (*F, n, Method::Auto());
        ok &= testNilpotentMinpoly (*F, n, Method::Elimination());
        ok &= testNilpotentMinpoly (*F, n, Method::Blackbox());
        typedef typename SparseMatrix<Field>::Row SparseVector;
        typedef DenseVector<Field> DenseVector;
        RandomDenseStream<Field, DenseVector, typename Field::NonZeroRandIter> zv_stream (*F, NzG, n, numVectors);
        RandomSparseStream<Field, SparseVector, typename Field::NonZeroRandIter > zA_stream (*F, NzG, (double) k / (double) n, n, n);

        ok &= testRandomMinpoly    (*F, iter, zA_stream, zv_stream, Method::Auto());
        ok &= testRandomMinpoly    (*F, iter, zA_stream, zv_stream, Method::Elimination());
        ok &= testRandomMinpoly    (*F, iter, zA_stream, zv_stream, Method::Blackbox());
        if (card>0){
            ok &= testGramMinpoly      (*F, n, Method::Auto());
            ok &= testGramMinpoly      (*F, n, Method::Elimination());
            ok &= testGramMinpoly      (*F, n, Method::Blackbox());
        }
            /*
              if(!ok)
              cout<<"FAILED"<<endl;
              else
              cout<<"PASS"<<endl;
            */

		delete F;
		nbiter--;
	}

	return ok;
}

int main (int argc, char **argv)
{
	commentator().setMaxDetailLevel (-1);

	commentator().setMaxDepth (-1);
	bool pass = true;

	integer q = -1;
	size_t b = 0; // set to a non zero value to force the bitsize of q
    int e = 1; // exponent for non prime fields
	size_t n = 30;
	int iterations = 3;
	int numVectors = 1;
	int k = 3;
    uint64_t seed = time(NULL);

	static Argument args[] = {
        { 'q', "-q Q", "Operate over the \"field\" GF(Q^E) [1].", TYPE_INTEGER, &q },
        { 'b', "-b B", "Set the bitsize of the field characteristic.", TYPE_INT, &b },
		{ 'e', "-e E", "Operate over the \"field\" GF(Q^E) [1].", TYPE_INT, &e },
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 'v', "-v V", "Use V test vectors for the random minpoly tests.", TYPE_INT,     &numVectors },
		{ 'k', "-k K", "K nonzero Elements per row in sparse random apply test.", TYPE_INT,     &k },
		{ 's', "-s seed", "set seed for the random generator.", TYPE_INT, &seed },
		END_OF_ARGUMENTS
	};


	parseArguments (argc, argv, args);
//	commentator().getMessageClass (TIMING_MEASURE).setMaxDepth (10);
//	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
//	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

    commentator().start ("Testing suite for minpoly", "testMinpoly", 1);

    pass &= run_with_field<Givaro::Modular<double> >(q,e,b,n,iterations,numVectors,k,seed);
    pass &= run_with_field<Givaro::Modular<int32_t> >(q,e,b,n,iterations,numVectors,k,seed);
    pass &= run_with_field<Givaro::Modular<Givaro::Integer> >(q,e,b?b:128,n/3+1,iterations,numVectors,k,seed);
    //pass &= run_with_field<Givaro::GFqDom<int64_t> >(q,e,b,n,iterations,numVectors,k,seed);
    pass &= run_with_field<Givaro::ZRing<Givaro::Integer> >(0,e,b?b:128,n/3+1,iterations,numVectors,k,seed);

    commentator().stop(MSG_STATUS(pass),(const char *) 0,"testMinpoly");

    return !pass;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
