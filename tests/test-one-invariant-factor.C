/* Copyright (C) LinBox
 *
 *  bds evolved this from test-last-invariant-factor
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



/*! @file  tests/test-one-invariant-factor.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */



#include <linbox/linbox-config.h>
#include "givaro/zring.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/ring/modular.h"
#include "linbox/algorithms/last-invariant-factor.h"
#include "linbox/algorithms/one-invariant-factor.h"
#include "linbox/blackbox/scompose.h"
#include "linbox/blackbox/random-matrix.h"
#include "linbox/algorithms/rational-solver.h"
#include <time.h>

#include "linbox/util/commentator.h"
#include "linbox/vector/stream.h"
#include "test-common.h"

using namespace LinBox;

template <class Ring, class OIF, class Vector>
bool testRandom(const Ring& R,
		OIF& oif,
		Vector& d)
		//LinBox::VectorStream<Vector>& stream1)
{

	std::ostringstream str;

	str << "Testing one invariant factor:";

        //commentator().start (str.str ().c_str (), "testRandom", stream1.m ());
        commentator().start (str.str ().c_str (), "testRandom", d.size());

        bool ret = true;

        VectorDomain<Ring> VD (R);

//	Vector d(R);

	typename Ring::Element x;

//	VectorWrapper::ensureDim (d, stream1.n ());

	size_t n = d. size();

//	 while (stream1) {

//		 commentator().startIteration ((unsigned)stream1.j ());

		 std::ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

 //               bool iter_passed = true;

//		stream1.next (d);

		report << "Input vector:  ";
		VD.write (report, d);
                report << endl;

		BlasMatrix<Ring> D(R, n, n), L(R, n, n), U(R, n, n), A(R,n,n),
		L2(R, n, n), U2(R, n, n);

		size_t i, j;

		for(i = 0; i < n; ++i) {
			R. assign (D[(size_t)i][(size_t)i], d[(size_t)i]);
			R. assign (L[(size_t)i][(size_t)i], R.one);
			R. assign (L2[(size_t)i][(size_t)i], R.one);
			R. assign (U[(size_t)i][(size_t)i], R.one);
			R. assign (U2[(size_t)i][(size_t)i], R.one);}

		for (i = 0; i < n; ++ i)

			for (j = 0; j < i; ++ j) {

				R.init(L[(size_t)i][(size_t)j], int64_t(rand() % 10));
				R.init(L2[(size_t)i][(size_t)j], int64_t(rand() % 10));

				R.init(U[(size_t)j][(size_t)i], int64_t(rand() % 10));
				R.init(U2[(size_t)j][(size_t)i], int64_t(rand() % 10));
			}


		BlasVector<Ring> tmp1(R,(size_t)n), tmp2(R,(size_t)n), e(R,(size_t)n);

		typename BlasMatrix<Ring>::ColIterator col_p;

		i = 0;
		for (col_p = A.colBegin();
		     col_p != A.colEnd(); ++ col_p, ++ i) {

			R.assign(e[(size_t)i],R.one);
			L2.apply(tmp2, e);
			U.apply(tmp1, tmp2);
			D.apply(tmp2, tmp1);
			// LinBox::BlasSubvector<BlasVector<Ring> > col_p_v(R,*col_p);
			// L.apply(col_p_v, tmp2);
			L.apply(tmp1, tmp2);
			U2.apply(*col_p, tmp1);
			R.assign(e[(size_t)i],R.zero);
		}


	for (size_t k = 0; k < n; ++k) {
		//size_t k = n/2;

		oif. oneInvariantFactor (x, A, k+1);


		report << "Computed one invariant factor: ";

		R. write (report, x);

		report << '\n';


		// typename BlasVector<Ring>::iterator p1;

		typename Ring::Element l;

		R. assign(l , d[k]);

//		for (p1 = d.begin(); p1 != d.end(); ++ p1)

//			R. lcmin (l, *p1);


		report << "Expected one invariant factor: ";


		R. write (report, l);

		report << '\n';

		if (!R. areEqual (l, x))

			ret = false;
	}

                if (!ret)

                        commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed one invariant factor is incorrect" << endl;



     //           commentator().stop ("done");

	 //}

	 //stream1.reset ();

	  commentator().stop (MSG_STATUS (ret), (const char *) 0, "testRandom");

	  return ret;

}

int main(int argc, char** argv)
{


        bool pass = true;

        static size_t n = 9;

	static unsigned int iterations = 1;

        static Argument args[] = {
                { 'n', "-n N", "Set order of test matrices to N.", TYPE_INT,     &n },
                { 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		END_OF_ARGUMENTS

        };

	parseArguments (argc, argv, args);

        typedef Givaro::ZRing<Integer>      Ring;

        Ring R; Ring::RandIter gen(R);

	commentator().start("One invariant factor test suite", "OIF");

        commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);


	typedef DixonSolver<Ring, Givaro::Modular<int32_t>, PrimeIterator<IteratorCategories::HeuristicTag> > Solver;
        // typedef DixonSolver<Ring, Givaro::Modular<double>, LinBox::PrimeIterator<IteratorCategories::HeuristicTag> > Solver;

	typedef LastInvariantFactor<Ring, Solver> LIF;
	typedef OneInvariantFactor<Ring, LIF, SCompose, RandomMatrix>  OIF;


	OIF oif;

	oif.  setThreshold (30);

	// Note:  oneInvariantFactor (used in smith-form-binary) is not designed to catch small prime factors
	vector<Ring::Element> d(n);
	d[0] = 1;
	for (size_t i = 1; i < n/2; ++i) d[i] = 101*d[i-1];
	for (size_t i = n/2; i < n; ++i) d[i] = d[i-1];
	d[n-3] *= 2*103;
	d[n-2] = 0;
	d[n-1] = 0;

    //RandomDenseStream<Ring> s1 (R, gen, n, iterations);
	if (!testRandom(R, oif, d)) pass = false;

    //RandomDenseStream<Ring> s2 (R, gen, n, iterations);
	//if (!testRandom(R, oif, s2, n/2)) pass = false;

	commentator().stop("One invariant factor test suite");
        return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
