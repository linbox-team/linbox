/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/solutions/det.h
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 * Modified by Pascal Giorgi and Zhendong Wan 2005/03
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

#include <linbox/integer.h>
#include <linbox/solutions/det.h>
#include <linbox/algorithms/cra.h>
#include <linbox/algorithms/matrix-mod.h>
#include <linbox/field/modular-double.h>
#include <linbox/randiter/random-prime.h>


namespace LinBox {


	/** Compute the determinant over {\bf Z} 
	 *
	 * Compute the determinant of a matrix.
	 * Perform the necessary modular reductions and
	 * reconstruct the result via Chinese remaindering
	 *
	 * @param res Element into which to store the result
	 * @param factory \ref{BlacboxFactory} that represents the matrix
	 */


	template <class Blackbox, class MethodTraits>
	typename Blackbox::Field::Element  &cra_det (typename Blackbox::Field::Element           &res,
						     const Blackbox                               &A,
						     const MethodTraits                           &M)
	{
		linbox_check (A.rowdim () == A.coldim ());
		
		commentator.start ("Determinant", "det");
		typedef Modular<double> Field;

		const unsigned int NUM_BITS = 21;
		integer start = 1 << (NUM_BITS);
				
		//int num_primes;
		//num_primes = ((int) (B.bitsize() + 2) /  NUM_BITS )  + (((int) (B.bitsize() + 2) % NUM_BITS)&1);

		//commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
		//	<< "Number of primes required: " << num_primes << endl;
		//cout<< "Number of primes required: " << num_primes << endl;

		//commentator.progress (0, num_primes);

		// I'm constructing all the fields at once in anticipation of
		// some parallelization infrastructure being implemented in the
		// not-too-distant future
		
				
		integer tmp;
		
		typedef typename MatrixModTrait<Blackbox, Field>::value_type FBlackbox;
		typename Field::Element  d;
		integer p;
		
		RandomPrime genprime(NUM_BITS);
		CRA<integer> cra;
		cra.initialize(1,1);
		
		while (!cra.terminated()) {
		
			genprime.randomPrime(p);			
			
			Field F(p);
			FBlackbox* Ap;
			MatrixMod::mod(Ap,A,F);

			det (d, *Ap, M);
			tmp=d;

			delete Ap;

			commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
				<< "Determinant modulo " << p << " is " << tmp << endl;
			commentator.progress ();

			cra.progress(p, tmp);
			
		}
		commentator.stop ("done", NULL, "det");
		
	
		cra.result(tmp);
		A.field().init(res,tmp);
		
		return res;
	}
}
