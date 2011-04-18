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


//#include <linbox/field/ntl-ZZ.h>
#include <linbox/field/PID-integer.h>
#include <time.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/util/commentator.h>
#include <linbox/vector/stream.h>
#include "test-common.h"
#include <linbox/blackbox/dense.h>
#include <linbox/algorithms/smith-form-adaptive.h>


template <class Ring, class SmithForm, class Vector>
bool testRandom(const Ring& R, 
		const SmithForm& SF,
		LinBox::VectorStream<Vector>& stream1) 
{
 
	
	std::ostringstream str;
        
	str << "Testing the adaptive algorithm for Smith form computation:\n";

        commentator.start (str.str ().c_str (), "testRandom");//, stream1.m ());

        bool ret = true;
        bool iter_passed = true;

        VectorDomain<Ring> VD (R);

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

		
		
		std::vector<integer> xi(A. rowdim());

		SF.smithForm (xi, A);
		typename Vector::iterator x_p; std::vector<integer>::iterator xi_p;
		for (x_p = x. begin(), xi_p = xi. begin(); x_p != x. end(); ++ x_p, ++ xi_p)
			A. field (). init (*x_p, *xi_p);
       
		
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
	  	  
	  commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandom");
                                                                                                        
	  return ret;

}

int main(int argc, char** argv) 
{
                                                                                                        
	bool pass = true;
	static size_t n = 35; 
	static int iterations = 1;
	static Argument args[] = {
		{ 'n', "-n N", "Set order of test matrices to N.", TYPE_INT,  &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		{ '\0' }
		};

	parseArguments (argc, argv, args);
	SmithFormAdaptive sf;

	commentator.start("Smith form adaptive algorithm test suite", "EGV++");
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);

	//typedef NTL_ZZ Ring; Ring R;
	typedef PID_integer Ring; Ring R;
	RandomDenseStream<Ring> s1 (R, n, iterations);
	pass = testRandom(R, sf, s1); 

	typedef PID_integer Ring2; Ring2 S;
	RandomDenseStream<Ring2> s2 (S, n, iterations);
	pass = pass && testRandom(S, sf, s2); 

	commentator.stop(MSG_STATUS(pass));
	return pass ? 0 : -1;
                                                                                                        
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
