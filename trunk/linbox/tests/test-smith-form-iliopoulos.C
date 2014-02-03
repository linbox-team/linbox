/* Copyright (C) LinBox
 *
 *  Author: Zhendong Wan
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


/*! @file   tests/test-smith-form-iliopoulos.C
 * @ingroup tests
 * @brief no doc.
 * @test no doc.
 */



//#include "linbox/field/ntl.h"
//#include "linbox/field/modular.h"
#include "linbox/field/PID-integer.h"
#include "linbox/field/PIR-ntl-ZZ_p.h"
#include "linbox/field/PIR-modular-int32.h"
//#include "linbox/integer.h"
#include "linbox/randiter/random-prime.h"
//#include "linbox/algorithms/last-invariant-factor.h"
#include "linbox/algorithms/smith-form-iliopoulos.h"
#include "linbox/algorithms/blas-domain.h"
//#include "linbox/algorithms/rational-solver.h"
//#include <time.h>
#include "linbox/util/commentator.h"
//#include "linbox/vector/stream.h"
#include "test-common.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/solutions/det.h"
#include <iostream>

#ifndef __LINBOX_HAVE_NTL
#error "you can't compile this test without NTL enabled. Please make sure you configured Linbox with --with-ntl=path/to/ntl"
#endif



using namespace LinBox;

template <class Ring>
bool testRead(const Ring& R, string file) {
	BlasMatrix<Ring> A(R);
	BlasMatrixDomain<Ring> BMD(R);
	std::ifstream data(file);
	A.read(data);
	size_t m = A.rowdim();
	size_t n = A.coldim();
	BlasMatrix<Ring> U(R, m, m), V(R, n, n), B(R, m, n);
	/* random unimodular
	BMD.randomNonsingular(U);
	BMD.randomNonsingular(V);
	*/
	for (size_t i = 0; i < m; ++i)
	{	for (size_t j = 0; j < m; ++j)
			U.setEntry(i, j, R.zero);
		U.setEntry(i, i, R.one);
	}
	for (size_t i = 0; i < n; ++i)
	{	for (size_t j = 0; j < n; ++j)
			V.setEntry(i, j, R.zero);
		V.setEntry(i, i, R.one);
	}
	BMD.copy(B, A);
	BMD.mulin_left(B,V);
	BMD.mulin_right(U,B);
	SmithFormIliopoulos::smithFormIn (B);
	return BMD.areEqual(A, B);
}

template <class Ring>
bool testRandom(const Ring& R, size_t n)
{
        bool pass = true;

        commentator().start ("Testing Iloipoulos elimination:", "testRandom");
	using namespace std;
	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	BlasMatrixDomain<Ring> BMD(R);
	BlasMatrix<Ring> D(R, n, n), L(R, n, n), U(R, n, n), A(R,n,n);

	// D is Smith Form
	const int m = 16;
	int p[m] = {1,1,1,1,1,1,2, 2, 2, 4, 3, 3, 3, 5, 7, 101};
	typename Ring::Element x, y;
	R.assign(x, R.one);
	if (n > 0) D.setEntry(0,0,x);
	for(size_t i = 1; i < n; ++i){
		R.init(y, p[rand()%m]);
		D.setEntry(i,i, R.mulin(x, y));
	}
	if (n > 0) D.setEntry(n-1,n-1, R.zero);

	// L and U are random unit triangular.
	for (size_t i = 0; i < n; ++ i) {
		for (size_t j = 0; j < i; ++ j) {
			L.setEntry(i,j, R.init(x, rand() % 10));
			U.setEntry(j,i, R.init(x, rand() % 10));
		}
		L.setEntry(i,i, R.one);
		U.setEntry(i,i, R.one);
	}

	// A is ULDU
	BMD.mul(A, U, L);
	BMD.mulin_left(A, D);
	BMD.mulin_left(A, U);

	D.write( report << "Smith form matrix:\n  " ) << endl;
	A.write( report << "input matrix:\n " ) << endl;
//	SmithFormIliopoulos::smithFormIn (A);
//	A.write( report << "smith of input matrix direct:\n " ) << endl;

	report << "Using PIRModular<int32_t>\n";
	typename Ring::Element d; R.init(d,16*101);//16*101*5*7*9);
	for (int i = 0; i < n-1; ++i) R.mulin(d, D.getEntry(x, i, i));
	R.write(report << "modulus: ", d) << endl;
	//det(d, D);
	//PIRModular<int32_t> Rd( (int32_t)(s % LINBOX_MAX_MODULUS));
	PIRModular<int32_t> Rp( (int32_t)d);
	BlasMatrix<PIRModular<int32_t> > Ap(Rp, n, n), Dp(Rp, n, n);
	MatrixHom::map (Ap, A);

	SmithFormIliopoulos::smithFormIn (Ap);

	Ap.write( report << "Computed Smith form: \n") << endl;
	MatrixHom::map (Dp, D);
	BlasMatrixDomain<PIRModular<int32_t> > BMDp(Rp);
	pass = pass and BMDp.areEqual(Dp, Ap);

	commentator().stop (MSG_STATUS (pass), (const char *) 0, "testRandom");
	return pass;

}

int main(int argc, char** argv)
{

        using namespace LinBox;

        bool pass = true;
        static size_t n = 2;
        static Argument args[] = {
                { 'n', "-n N", "Set order of test matrices to N.", TYPE_INT,     &n },
		END_OF_ARGUMENTS
        };
        parseArguments (argc, argv, args);

	commentator().start("Ilioloulos Smith Form test suite", "Ilioloulos");
        commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);

        PID_integer R; // Ring of integers
        //NTL_ZZ R;

        if (!testRandom(R, n)) pass = false;
//        if (!testRead(R, "data/Ismith.mat")) pass = false;

	commentator().stop("Ilioloulos Smith Form test suite");
        return pass ? 0 : -1;

}


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
