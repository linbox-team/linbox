/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/solutions/det.h
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __DET_H
#define __DET_H


#include "linbox/field/modular.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/blackbox/factory.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h"
#include "linbox/algorithms/cra.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/solutions/methods.h"
#include "linbox/util/prime-stream.h"
#include "linbox/util/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
	/** Determinant over a field
	 *
	 * Compute the determinant of a linear operator A, represented as a
	 * black box, over a field F.
	 *
	 * This implementation is essentially direct, in that it does not
	 * perform any modular reduction and reconstruction. Thus, it is not
	 * recommended that one use this function to compute the determinant of
	 * an integer or rational matrix. One should instead use the version
	 * indicated below that uses \ref{BlackboxFactory}.
	 *
	 * @param res Field element into which to store the result
	 * @param A Black box of which to compute the determinant
	 * @param F Field over which to compute the determinant
	 * @param M Method traits
	 */

	template <class Field, class Blackbox>
	typename Field::Element &det (typename Field::Element         &res,
				      const Blackbox                  &A,
				      const Field                     &F,
				      const MethodTrait::Wiedemann    &M = MethodTrait::Wiedemann ()) 
	{
		typedef std::vector<typename Field::Element> Polynomial;

		commentator.start ("Determinant", "det");

		linbox_check (A.coldim () == A.rowdim ());

		Polynomial               phi;
		unsigned long            deg;
		typename Field::RandIter iter (F);

		// Precondition here to separate the eigenvalues, so that
		// minpoly (B) = charpoly (B) with high probability

		std::vector<typename Field::Element> d (A.coldim ());
		typename Field::Element pi;
		size_t i;
		size_t iternum = 1;
		do {

			F.init (pi, 1);

			for (i = 0; i < A.coldim (); i++) {
				do iter.random (d[i]); while (F.isZero (d[i]));
				F.mulin (pi, d[i]);
			}
	
			Diagonal<Field> D (F, d);

			Compose<Blackbox,Diagonal<Field> > B (&A, &D);

			typedef Compose<Blackbox,Diagonal<Field> > Blackbox1;

			BlackboxContainer<Field, Blackbox1> TF (&B, F, iter);
			
			MasseyDomain<Field, BlackboxContainer<Field, Blackbox1> > WD (&TF, M.earlyTermThreshold ());

			WD.minpoly (phi, deg);
			//cout << "\tdet: iteration # " << iternum << "\tMinpoly deg= " << phi.size() << "\n";
			
			++iternum;
		} while (!F.isZero (phi[0]) && phi.size () < A.coldim () + 1);

		if (deg & 1 == 1)
			F.negin (pi);

		F.div (res, phi[0], pi);

		commentator.stop ("done", NULL, "det");

		return res;
	}

// 	template <class Field, class Blackbox>
// 	typename Field::Element &det (typename Field::Element         &res,
// 				      const Blackbox                  &A,
// 				      const Field                     &F,
// 				      const MethodTrait::Wiedemann    &M = MethodTrait::Wiedemann ()) 
// 	{
// 		return det<Field, Blackbox, std::vector<typename Field::Element> >(res,A,F,M);
// 	}
	

	/** Determinant over $\mathbb{Z}$ or $\mathbb{Q}$
	 *
	 * Compute the determinant of a matrix, represented via a
	 * \ref{BlackboxFactory}. Perform the necessary modular reductions and
	 * reconstruct the result via Chinese remaindering or rational number
	 * reconstruction.
	 *
	 * @param res Element into which to store the result
	 * @param factory \ref{BlacboxFactory} that represents the matrix
	 */

	// FIXME: Right now we only support doing this over Modular<uint32> --
	// that's probably a bad idea. There needs to be a way to get from the
	// field some idea of where a "good" choice of moduli is.

	template <class Element, class Blackbox, class Field>
	Element &det (Element                &res,
		      BlackboxFactory<Field, Blackbox> &factory) 
	{
		linbox_check (factory.rowdim () == factory.coldim ());

		commentator.start ("Determinant", "det");

		// Number of bits -- note that this assumes we are in Modular<uint32>
		const unsigned int NUM_BITS = 31;

		integer start = 1 << (NUM_BITS - 1);
		PrimeStream<typename Field::Element> stream (start);

		// Get the log in base 2^(NUM_BITS - 1) of the Hadamard bound on the
		// determinant of the matrix
		integer B;
		double n = factory.rowdim ();
		int num_primes;

		factory.maxNorm (B);

		// If this overflows an integer, the problem is just impossible
		// anyway, so I'm assuming it won't.
		num_primes = (int) ceil (n * (log (n) / 2.0 + log ((double) B)) / (M_LN2 * (double) (NUM_BITS - 1)));

		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "Number of primes required: " << num_primes << endl;

		commentator.progress (0, num_primes);

		// I'm constructing all the fields at once in anticipation of
		// some parallelization infrastructure being implemented in the
		// not-too-distant future
		typename std::vector<Field> F;
		std::vector<uint32> moduli (num_primes);
		std::vector<uint32> res_mod (num_primes);

		for (int i = 0; i < num_primes; ++i) {
			stream >> moduli[i];
			F.push_back (Field (moduli[i]));
		}

		// In principal, the following could be done in parallel. I'm
		// isolating the loop to make that easier. I envision something
		// like this:
		//
		// Field::Element &cb (Field::Element &res,
		//                     Field &F,
		//                     BlackboxFactory<Field> &factory)
		// {
		//         ... do individual computation ...
		//         return res;
		// }
		//
		// Element &det (...)
		// {
		//         ...
		//         ThreadManager manager;
		//         manager.run (cb, data, results);
		// }
		//
		// The object ThreadManager is responsible for (1) knowing about
		// the parallelization capabilities of the machine, and (2)
		// invoking the callback on the optimal number of threads given
		// the architecture. It also runs the callback in serial as
		// necessary.
		//
		// Anyway, back to coding...

		for (int i = 0; i < num_primes; ++i) {
			Blackbox *A = factory.makeBlackbox (F[i]);
			det (res_mod[i], *A, F[i]);
			delete A;

			commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
				<< "Determinant modulo " << moduli[i] << " is " << res_mod[i] << endl;

			commentator.progress ();
		}

		cra (res, res_mod, moduli);

		commentator.stop ("done", NULL, "det");
		return res;
	}
}

#endif // __DET_H
