/* Copyright (C) LinBox
 *
 * Author: Zhendong Wan
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



#include <linbox/field/PID-integer.h>
#include <linbox/field/modular-int32.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/diagonal.h>
#include <linbox/algorithms/rational-solver.h>
#include <linbox/randiter/random-prime.h>
#include <iostream>
#include "test-common.h"
#include "linbox/vector/stream.h"
#include "linbox/util/commentator.h"

/// Testing Nonsingular Random Diagonal solve.
template <class Ring, class Field, class Vector>
bool testRandomSolve (const Ring& R,
		      const Field& f,
		      LinBox::VectorStream<Vector>& stream1,
		      LinBox::VectorStream<Vector>& stream2) 
{

	std::ostringstream str;
	
	commentator.start ("Testing Nonsingular Random Diagonal solve ","testNonsingularRandomDiagonalSolve", stream1.size());// "testNonsingularRandomMatrixSolve", stream1.m ());

	bool ret = true;

        bool iter_passed = true;
	
	VectorDomain<Ring> VD (R);

	Vector d, b, x, y;

	VectorWrapper::ensureDim (d, stream1.n ());
        VectorWrapper::ensureDim (b, stream1.n ());
        VectorWrapper::ensureDim (x, stream1.n ());
        VectorWrapper::ensureDim (y, stream1.n ());

	int n = d. size();

	while (stream1 && stream2) {
        
		commentator.startIteration (stream1.j ());
                                                                                                        
                //ActivityState state = commentator.saveActivityState ();
                                                                                                        
                iter_passed = true;
                
		bool zeroEntry;
		do {
		  stream1.next (d);
		  zeroEntry = false;
		  for (size_t i=0; i<stream1.n(); i++)
		    zeroEntry |= R.isZero(d[i]);
		} while (zeroEntry);
		
                stream2.next (b);
                                                                                                        
                std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
                report << "Diagonal entries: ";
                VD.write (report, d);
                report << endl;
                                                                                                        
                report << "Right-hand side:  ";
                VD.write (report, b);
                report << endl;

                //Diagonal<Ring> D(R, d);                                                                               
		
		DenseMatrix<Ring> D(R, n, n);
	 
		for(int i = 0; i < n; ++i) R.init (D[i][i],  d[i]);
						
		typedef RationalSolver<Ring, Field, LinBox::RandomPrimeIterator> RSolver;
		RSolver rsolver;
 
		//std::vector<std::pair<typename Ring::Element, typename Ring::Element> > answer(n);
		std::vector<typename Ring::Element> num(n);
		typename Ring::Element den;
 
		SolverReturnStatus solveResult = rsolver.solve(num, den, D, b, 30); //often 5 primes are not enough
		
		/*
		typename Ring::Element lden;

		R. init (lden, 1);

		typename std::vector<std::pair<typename Ring::Element, typename Ring::Element> >::iterator p;
		
		for (p = answer.begin(); p != answer.end(); ++ p)
			R. lcm (lden, lden, p->second);
		typename Vector::iterator p_x;
		//typename Vector::iterator p_y;
		*/

		if (solveResult == SS_OK) {
		/*
		  for (p = answer.begin(), p_x = x. begin(); 
		       p != answer.end();
		       ++ p, ++ p_x) {
		    
		    R. mul (*p_x, p->first, lden);
		    
		    R. divin (*p_x, p->second);
		    
		  }
		  
		  D. apply (y, x);
		  */
		  D. apply (y, num);
		  
		  VD. mulin(b, den);
		  
		  if (!VD.areEqual (y, b)) {
		    ret = iter_passed = false;
		    commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
		      << "ERROR: Computed solution is incorrect" << endl;
		  }
		}
		else {
		    ret = iter_passed = false;
		    commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
		      << "ERROR: Did not return OK solving status" << endl;		  
		}
		
		commentator.stop ("done");
                commentator.progress ();
		
	}

	
	stream1.reset ();
        stream2.reset ();
	
        commentator.stop (MSG_STATUS (ret), (const char *) 0, "testNonsingularRandomDiagonalSolve");

	return ret;
}	

int main(int argc, char** argv) 
{

	
	bool pass = true;
 
        static size_t n = 10;
                
	static int iterations = 1;
 
        static Argument args[] = {
                { 'n', "-n N", "Set order of test matrices to N.", TYPE_INT, &n},
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		{ '\0' }
        };


	parseArguments (argc, argv, args);
	
	typedef Modular<LinBox::int32> Field;
	
	typedef PID_integer     Ring;

	Ring R;

	Field F(101);
	
	RandomDenseStream<Ring> s1 (R, n, iterations), s2 (R, n, iterations);

	if (!testRandomSolve(R, F, s1, s2)) pass = false;
	
	return pass ? 0 : -1;
	
}

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
