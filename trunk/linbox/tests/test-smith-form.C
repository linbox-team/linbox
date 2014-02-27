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

/*! @file tests/test-smith-form.C
 * @ingroup tests
 * @brief no doc.
 * @test no doc.
	//!@bug should work for NTL integers too
 */


#include <time.h>
//#ifdef __LINBOX_HAVE_NTL
//#include "linbox/field/ntl-ZZ.h"
//#endif
#include "linbox/field/PID-integer.h"
#include "linbox/util/commentator.h"
#include "linbox/vector/stream.h"
#include "test-common.h"
#include "linbox/vector/blas-vector.h"
#include "linbox/solutions/smith-form.h"
using LinBox::parseArguments;
using LinBox::commentator;
using LinBox::Commentator;
using LinBox::integer;
using LinBox::PID_integer;
using LinBox::BlasMatrix;
using LinBox::BlasVector;

template <class Ring, class Vector>
bool testRandom(const Ring& R,
		LinBox::VectorStream<Vector>& stream1)
{

	std::ostringstream str;

	str << "Testing the smithForm function in solutions directory:\n";

        commentator().start (str.str ().c_str (), "testRandom");//, stream1.m ());

        bool ret = true;

        LinBox::VectorDomain<Ring> VD (R);

	Vector d(R), x(R);

	LinBox::VectorWrapper::ensureDim (d, stream1.n ());

	LinBox::VectorWrapper::ensureDim (x, stream1.n ());


	int n = (int)d. size();

	 while (stream1) {

                commentator().startIteration ((unsigned)stream1.j ());

		std::ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

                bool iter_passed = true;

                stream1.next (d);

		report << "Input vector:  ";
		VD.write (report, d);
                report << endl;

		BlasMatrix<Ring> D(R, n, n), L(R, n, n), U(R, n, n), A(R,n,n);

		int i, j;

		for(i = 0; i < n; ++i) {
			R. assign (D[(size_t)i][(size_t)i], d[(size_t)i]);
			R. assign(L[(size_t)i][(size_t)i], R.one);
			R. assign(U[(size_t)i][(size_t)i], R.one);}

		for (i = 0; i < n; ++ i)

			for (j = 0; j < i; ++ j) {

				R.init(L[(size_t)i][(size_t)j], rand() % 10);

				R.init(U[(size_t)j][(size_t)i], rand() % 10);
			}


		BlasVector<Ring> tmp1(R,(size_t)n), tmp2(R,(size_t)n), e(R,(size_t)n);

		typename BlasMatrix<Ring>::ColIterator col_p;

		i = 0;
		for (col_p = A.colBegin();
		     col_p != A.colEnd(); ++ col_p, ++ i) {

			R.assign(e[(size_t)i],R.one);
			U.apply(tmp1, e);
			D.apply(tmp2, tmp1);
			// LinBox::BlasSubvector<BlasVector<Ring> > col_p_v (R, *col_p);
			// L.apply(col_p_v, tmp2);
			L.apply(*col_p, tmp2); //! @internal @bug  should use Triangular apply ? We are doing this many times, factor somewhere in test-utils.h ? why not some ftrtr routine for that ?
			R.assign(e[(size_t)i],R.zero);
		}

		typename Vector::iterator x_p;
		PID_integer Z;
		BlasVector<PID_integer> xi(Z,A. rowdim());
		BlasVector<PID_integer>::iterator xi_p;
		std::list<std::pair<integer, size_t> > cpt;
		smithForm (cpt, A);
		std::list<std::pair<integer, size_t> >::iterator cpt_p;

		xi_p = xi. begin();
		for (cpt_p = cpt.begin(); cpt_p != cpt.end(); ++ cpt_p) {
			for (size_t ii = 0; ii < cpt_p -> second; ++ ii, ++ xi_p)
				*xi_p = cpt_p -> first;
		}

		for (x_p = x. begin(), xi_p = xi. begin(); x_p != x. end(); ++ x_p, ++ xi_p)
			A. field (). init (*x_p, *xi_p);

		report << "Computed Smith form: \n";
		VD. write (report, x);

		report << '\n';

		typename BlasVector<Ring>::iterator p1, p2;
		typename Ring::Element g;

		for (p1 = d.begin(); p1 != d.end(); ++ p1) {
			for ( p2 = p1 + 1; p2 != d.end(); ++ p2) {
				if (R. isUnit(*p1))  break;
				else if (R. isZero (*p2)) continue;
				else if (R. isZero (*p1)) std::swap (*p1, *p2);
				else { // (*p1, *p2) <-- (g, *p1 * *p2 / g), where g = gcd(*p1, *p2)
					R. gcd (g, *p1, *p2);
					R. divin (*p2, g);
					R. mulin (*p2, *p1);
					R. assign (*p1, g);
				}
			}
		}

		report << "Expected smith form:\n";
		VD.write (report, d) << endl;

		if (!VD.areEqual (d, x))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed Smith form is incorrect" << endl;

		commentator().stop ("done");
		commentator().progress ();
	 }

	//stream1.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testRandom");

	return ret;
}

int main(int argc, char** argv)
{

	bool pass = true;
	static size_t n =3;
	static int iterations = 2;
	static Argument args[] = {
		{ 'n', "-n N", "Set order of test matrices to N.", TYPE_INT,  &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);
	//!@bug should be tried on NTZ_LL too
	typedef LinBox::PID_integer      Ring;

	Ring R;

	commentator().start("Smith form test suite", "Smith");
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);

	LinBox::RandomDenseStream<Ring> s1 (R, n, (unsigned int)iterations);
	if (!testRandom(R, s1)) pass = false;

#if 0
#ifdef __LINBOX_HAVE_NTL
	typedef LinBox::NTL_ZZ      Ring2;
	Ring2 R2;

	LinBox::RandomDenseStream<Ring2> s2 (R2, n, (unsigned int)iterations);
	if (!testRandom(R2, s2)) pass = false;

#endif
#endif

	commentator().stop("Smith form test suite");
	return pass ? 0 : -1;

}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

