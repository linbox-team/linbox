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





#include <linbox/field/ntl-ZZ.h>
#include <linbox/field/modular-int32.h>
#include <linbox/field/PIR-ntl-ZZ_p.h>
#include <linbox/field/PIR-modular-int32.h>
#include <linbox/integer.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/blackbox/dense.h>
#include <linbox/algorithms/last-invariant-factor.h>
#include <linbox/algorithms/smith-form-iliopoulos.h>
#include <linbox/algorithms/rational-solver.h>
#include <time.h>
#include <linbox/util/commentator.h>
#include <linbox/vector/stream.h>
#include "test-common.h"
#include <linbox/algorithms/matrix-hom.h>

#define int32 LinBox::int32

using namespace LinBox;

template <class Ring, class Vector>
bool testRandom(const Ring& R, 
		LinBox::VectorStream<Vector>& stream1) 
{
 
	using namespace std;
	
	ostringstream str;
        
	str << "Testing Iloipoulos elimination:";

        commentator.start (str.str ().c_str (), "testRandom", stream1.m ());

        bool ret = true;
                                                                                                        
        bool iter_passed = true;
                                                                                                        
        VectorDomain<Ring> VD (R);

	Vector d, x;
	
	VectorWrapper::ensureDim (d, stream1.n ());
	
	VectorWrapper::ensureDim (x, stream1.n ());
	
	int n = d. size();

	 while (stream1) {
                                                                                                        
                commentator.startIteration (stream1.j ());
                                                                                                        
		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);  

                iter_passed = true;
                                                                                                        
                stream1.next (d);

		report << "Input vector:  ";
		VD.write (report, d);
                report << endl;

		DenseMatrix<Ring> D(R, n, n), L(R, n, n), U(R, n, n), A(R,n,n);

		int i, j;

		for(i = 0; i < n; ++i) { R. assign (D[i][i], d[i]); R. init (L[i][i], 1); R. init (U[i][i], 1);}
		
		for (i = 0; i < n; ++ i) 
		
			for (j = 0; j < i; ++ j) {
				
				R.init(L[i][j], rand() % 10);
				
				R.init(U[j][i], rand() % 10);
			}
	
	
		std::vector<typename Ring::Element> tmp1(n), tmp2(n), e(n);
		
		typename DenseMatrix<Ring>::ColIterator col_p;

		i = 0;
		for (col_p = A.colBegin(); col_p != A.colEnd(); ++ col_p, ++ i) {
			
			R.init(e[i],1);
			U.apply(tmp1, e);
			D.apply(tmp2, tmp1);
			L.apply(*col_p, tmp2);
			R.init(e[i],0);
		}

		
		
		typename Ring::Element s;
		
		R. init (s, 1);

		typename Vector::iterator d_p;

		for(d_p = d.begin(); d_p!= d.end(); ++ d_p)
			
			R. lcm (s, s, *d_p);

		report << "Last Invariant factor: " ;
		
		R. write (report, s);

		report << '\n';
		
		
		if (s >= LINBOX_MAX_MODULUS) {

			report << "Using PIR_ntl_ZZ_p\n";
			
			PIR_ntl_ZZ_p PIR(s);
			
			DenseMatrix<PIR_ntl_ZZ_p> Ap(PIR, A.rowdim(), A.coldim());
			
			MatrixHom::map (Ap, A, PIR);
			
			SmithFormIliopoulos::smithFormIn (Ap);
			
			report << "Computed Smith form: \n";
			
			for ( unsigned int i = 0; i < A. rowdim(); ++ i)
				report << Ap[i][i] << " ";
			
			report << '\n';
			
			int i = 0;
			
			typename std::vector<typename Ring::Element>::iterator p1;
			
			
			for (p1 = x. begin(); p1 != x. end(); ++ p1, ++ i) {
				
				if (PIR.isZero(Ap[i][i])) 
					
					R.assign (*p1, s);
				
				else
					
					R.assign (*p1, NTL::rep(Ap[i][i]));
			}
		}

		else {

			report << "Using PIRModular<int32>\n";
		
			PIRModular<int32> PIR(s % LINBOX_MAX_MODULUS);
			
			DenseMatrix<PIRModular<int32> > Ap(PIR, A.rowdim(), A.coldim());
			
			MatrixHom::map (Ap, A, PIR);
			
			SmithFormIliopoulos::smithFormIn (Ap);
			
			
			report << "Computed Smith form: \n";
 
			for ( unsigned int i = 0; i < A. rowdim(); ++ i)
				report << Ap[i][i] << " ";
			
			report << '\n';
			
			
			typename std::vector<typename Ring::Element>::iterator p1;
			int i = 0;
			
			for (p1 = x. begin(); p1 != x. end(); ++ p1, ++ i) {
				
				if (PIR.isZero(Ap[i][i]))
					
					R.assign (*p1, s);
				
				else
 
					R.init (*p1, Ap[i][i]);
			}
		}
			
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
		
		for (p1 = d.begin(); p1 != d.end(); ++ p1)
			report << * p1 << " ";

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
                                                                                                        
        using namespace LinBox;
                                                                                                        
        bool pass = true;
                                                                                                        
        static size_t n = 20;
                                                                                                        
        static int iterations = 1;
                                                                                                        
        static Argument args[] = {
                { 'n', "-n N", "Set order of test matrices to N.", TYPE_INT,     &n },
                { 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
				{ '\0' }
        };
                                                                                                        
                                                                                                        
        parseArguments (argc, argv, args);
                                                                                                        
        typedef NTL_ZZ      Ring;
                                                                                                        
        Ring R;

	commentator.start("Ilioloulos Smith Form test suite", "Ilioloulos");

        commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);

        RandomDenseStream<Ring> s1 (R, n, iterations);
                                                                                                        
        if (!testRandom(R, s1)) pass = false;
                                                                                                        
	commentator.stop("Ilioloulos Smith Form test suite");
        return pass ? 0 : -1;
                                                                                                        
}

#undef int32
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
