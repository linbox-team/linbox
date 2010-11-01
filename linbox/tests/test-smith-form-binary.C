/* Copyright (C) LinBox
 *
 *
 *  Author: Zhendong Wan
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




#include <linbox/field/PID-integer.h>
#ifdef __LINBOX_HAVE_NTL
#include "linbox/field/ntl-ZZ.h"
#endif
#include <linbox/field/modular-int32.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/blackbox/dense.h>
#include <linbox/algorithms/matrix-rank.h>
#include <linbox/algorithms/last-invariant-factor.h>
#include <linbox/algorithms/one-invariant-factor.h>
#include <linbox/algorithms/smith-form-binary.h>
#include <linbox/blackbox/scompose.h>
#include <linbox/blackbox/random-matrix.h>
#include <linbox/algorithms/rational-solver.h>
#include <time.h>

#include <linbox/util/commentator.h>
#include <linbox/vector/stream.h>
#include "test-common.h"

template <class Ring, class SmithForm, class Vector>
bool testRandom(const Ring& R, 
		const SmithForm& SF,
		LinBox::VectorStream<Vector>& stream1) 
{

    std::ostringstream str;
        
	str << "Testing Smith Form binary(EGV++):";

        commentator.start (str.str ().c_str (), "testSmithform");//, stream1.m ());

        bool ret = true;
        bool iter_passed = true;

        LinBox::VectorDomain<Ring> VD (R);

	Vector d, x;
	
	VectorWrapper::ensureDim (d, stream1.n ());
	
	VectorWrapper::ensureDim (x, stream1.n ());

	
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

		
		
		SF.smithFormBinary (x, A);
       
		
		report << "Computed Smith form: \n";
		
		VD. write (report, x);
		
		report << '\n';
	
				
		typename std::vector<typename Ring::Element>::iterator p1, p2;
		typename Ring::Element g;
		
		
		for (p1 = d.begin(); p1 != d.end(); ++ p1) {
			
			for ( p2 = p1 + 1; p2 != d.end(); ++ p2) {
				
				if (R. isUnit(*p1))  break;
 
				else if (R. isZero (*p2)) continue;
				
				else if (R. isZero (*p1)) {
                                                std::swap (*p1, *p2);
				}
				
				else {
					R. gcd (g, *p1, *p2);
					
					R. divin (*p2, g);
					
					R. mulin (*p2, *p1);
					
					R. assign (*p1, g);
				}
			}
		}
		

		report << "Expected smith form:\n";
		
		VD.write (report, d);

		report << '\n';

		if (!VD.areEqual (d, x))
			
			ret = iter_passed = false;
		
                if (!iter_passed) 
			
                        commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed Smith form is incorrect" << endl;
			
		

                commentator.stop ("done");

                commentator.progress ();
		
	 }
	 
	 //stream1.reset ();
	  	  
	  commentator.stop (MSG_STATUS (ret), (const char *) 0, "testSmithform");
                                                                                                        
	  return ret;

}

int main(int argc, char** argv) 
{

	bool pass = true;

	static size_t n =5; 

	static int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set order of test matrices to N.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.start("SmithFormBinary test suite", "SmithFormBinary");
	std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	{
		typedef PID_integer Ring;

		Ring R;

		report << std::endl << "EGV++ algorithm test suite with LinBox/Givaro PID:\n";

		commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);

		RandomDenseStream<Ring> s1 (R, n, iterations);

		typedef Modular<LinBox::int32> Field;
		typedef RationalSolver<Ring, Field, LinBox::RandomPrimeIterator> Solver;
		typedef LastInvariantFactor<Ring, Solver> LIF;
		typedef OneInvariantFactor<Ring, LIF, SCompose, RandomMatrix>  OIF;
		typedef SmithFormBinary<Ring, OIF, MatrixRank<Ring, Field > > SF;

		SF sf;
		sf. setOIFThreshold (30);
		sf. setLIFThreshold (30);

		if (!testRandom(R, sf, s1)) pass = false;
	}
        
#if 0
//#ifdef __LINBOX_HAVE_NTL
// NTL_ZZ not working here
	{
		typedef NTL_ZZ Ring;

		Ring R;

		report << std::endl << "EGV++ algorithm test suite with NTL_ZZ :\n";
		commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);

		RandomDenseStream<Ring> s1 (R, n, iterations);

		typedef Modular<LinBox::int32> Field;
		typedef RationalSolver<Ring, Field, LinBox::RandomPrimeIterator> Solver;
		typedef LastInvariantFactor<Ring, Solver> LIF;
		typedef OneInvariantFactor<Ring, LIF, SCompose, RandomMatrix>  OIF;
		typedef SmithFormBinary<Ring, OIF, MatrixRank<Ring, Field > > SF;

		SF sf;
		sf. setOIFThreshold (30);
		sf. setLIFThreshold (30);

		if (!testRandom(R, sf, s1)) pass = false;
	}
#endif

	commentator.stop("SmithFormBinary test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
