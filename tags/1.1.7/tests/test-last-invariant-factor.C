/* Copyright (C) LinBox
 *
 *  Author: Zhendong Wan
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


/** \file test-last-invariant-factor.C
 */

#include <linbox/field/PID-integer.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/field/modular-int32.h>
#include <linbox/blackbox/dense.h>
#include <linbox/algorithms/matrix-rank.h>
#include <linbox/algorithms/last-invariant-factor.h>
#include <linbox/blackbox/scompose.h>
#include <linbox/blackbox/random-matrix.h>
#include <linbox/algorithms/rational-solver.h>
#include <time.h>


#include <linbox/util/commentator.h>
#include <linbox/vector/stream.h>
#include "test-common.h"

template <class Ring, class LIF, class Vector>
bool testRandom(const Ring& R, 
		const LIF& lif,
		LinBox::VectorStream<Vector>& stream1) 
{
 
	std::ostringstream str;
        
	str << "Testing last invariant factor:";

        commentator.start (str.str ().c_str (), "testRandom", stream1.m ());

        bool ret = true;
        bool iter_passed = true;

        VectorDomain<Ring> VD (R);

	Vector d;

	typename Ring::Element x;
	
	VectorWrapper::ensureDim (d, stream1.n ());
	
	int n = d. size();

	 while (stream1) {
		 
		 commentator.startIteration (stream1.j ());
                 
		 std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);  

                iter_passed = true;
		
		stream1.next (d);

		report << "Input vector:  ";
		VD.write (report, d);
                report << endl;

		DenseMatrix<Ring> D(R, n, n), L(R, n, n), U(R, n, n), A(R,n,n);

		int i, j;

		for(i = 0; i < n; ++i) { 
			R. assign (D[i][i], d[i]); 
			R. init (L[i][i], 1);
			R. init (U[i][i], 1);}
		
		for (i = 0; i < n; ++ i) 
		
			for (j = 0; j < i; ++ j) {
				
				R.init(L[i][j], rand() % 10);
				
				R.init(U[j][i], rand() % 10);
			}
	
	
		std::vector<typename Ring::Element> tmp1(n), tmp2(n), e(n);
		
		typename DenseMatrix<Ring>::ColIterator col_p;

		i = 0;
		for (col_p = A.colBegin(); 
		     col_p != A.colEnd(); ++ col_p, ++ i) {
			
			R.init(e[i],1);
			U.apply(tmp1, e);
			D.apply(tmp2, tmp1);
			L.apply(*col_p, tmp2);
			R.init(e[i],0);
		}

		
		
		lif. lastInvariantFactor (x, A);
       
		
		report << "Computed last invariant factor: \n";
		
		R. write (report, x);
		
		report << '\n';
	
				
		typename std::vector<typename Ring::Element>::iterator p1;
		
		typename Ring::Element l;
		
		R. init (l , 1);
		
		for (p1 = d.begin(); p1 != d.end(); ++ p1) 
			
			R. lcmin (l, *p1);

			

		report << "Expected last invariant factor: \n";
		
		R. write (report, l);

		report << '\n';

		if (!R. areEqual (l, x))
			
			ret = iter_passed = false;
		
                if (!iter_passed) 
			
                        commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed last invariant factor is incorrect" << endl;
			
		

                commentator.stop ("done");

                commentator.progress ();
		
	 }
	 
	 //stream1.reset ();
	  	  
	  commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandom");
                     
	  return ret;

}

int main(int argc, char** argv) 
{
        

        bool pass = true;
        
        static size_t n = 10;
        
	static int iterations = 1;
        
        static Argument args[] = {
                { 'n', "-n N", "Set order of test matrices to N.", TYPE_INT,     &n },
                { 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
				{ '\0' }

        };
	
	parseArguments (argc, argv, args);
        
        typedef PID_integer      Ring;
                                        
        Ring R;

	commentator.start("Last invariant factor test suite", "LIF");

        commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);

        RandomDenseStream<Ring> s1 (R, n, iterations);

	typedef RationalSolver<Ring, Modular<LinBox::int32>, LinBox::RandomPrimeIterator> Solver;

	typedef LastInvariantFactor<Ring, Solver> LIF;

	LIF lif;
	
	lif.  setThreshold (30);

	if (!testRandom(R, lif, s1)) pass = false;
                              
	commentator.stop("Last invariant factor test suite");
        return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
