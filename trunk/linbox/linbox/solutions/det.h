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

#include "linbox/blackbox/blas-blackbox.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h"
#include "linbox/algorithms/cra.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/solutions/methods.h"
#include "linbox/util/prime-stream.h"
#include "linbox/util/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
	// for specialization with respect to the DomainCategory
	template< class Blackbox, class MethodTraits, class DomainCategory>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &res, 
						const Blackbox                              &A,
						const DomainCategory                      &tag,
						const MethodTraits                          &M);

	/** Compute the determinant of A
	 *
	 * The determinant of a linear operator A, represented as a
	 * black box, is computed over the ring or field of A.
	 *
	 * @param res Field element into which to store the result
	 * @param A Black box of which to compute the determinant
	 * @param M Method traits
	 */
	template <class Blackbox, class MethodTraits>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &res, 
						const Blackbox                              &A,				
						const MethodTraits                           &M) 
	{
		return det(res, A, FieldTraits<typename Blackbox::Field>::categoryTag(), M);
	}

	// The det with default MethodTrait 
	template<class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &res, 
						const Blackbox                               &A)
	{
		return det(res, A, FieldTraits<typename Blackbox::Field>(), MethodTrait::BlasElimination());
	}

	// The entry domain must be a field (and ought to be finite).
	template <class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &res, 
						const Blackbox                              &A,
						const RingCategories::ModularTag          &tag,
						const MethodTrait::Wiedemann                &M) 
	{
		typedef typename Blackbox::Field Field;
		typedef std::vector<typename Field::Element> Polynomial;
		Field F = A.field();
		
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



	// ring should be a (finite) field.
	template <class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &res,
						const Blackbox                              &A,
						const RingCategories::ModularTag          &tag,
						const MethodTrait::BlasElimination           &M) 
	{
		typedef typename Blackbox::Field Field;
		Field F = A.field();
		
		commentator.start ("Determinant", "det");

		linbox_check (A.coldim () == A.rowdim ());

		BlasMatrix<typename Field::Element> B(A);
		BlasMatrixDomain<Field> BMD(F);
		res= BMD.det(B);
		commentator.stop ("done", NULL, "det");

		return res;
	}

	
	template <class Field>
	typename Field::Element &detin (typename Field::Element             &res,
					BlasBlackbox<Field>                   &A,
					const MethodTrait::BlasElimination     &M) 
	{
		Field F = A.field();
		
		commentator.start ("Determinant", "det");
		linbox_check (A.coldim () == A.rowdim ());

		BlasMatrixDomain<Field> BMD(F);
		res= BMD.detin(static_cast<BlasMatrix<typename Field::Element>& > (A));
		commentator.stop ("done", NULL, "det");

		return res;
	}
} // end of LinBox namespace 

#include "linbox/algorithms/cra-det-integer.h"

namespace LinBox {
	
	template <class Blackbox, class MethodTrait>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &res,
						const Blackbox                              &A,
						const RingCategories::IntegerTag          &tag,
						const MethodTrait                           &M) 
	{
		return cra_det (res, A, M);
	}


} // end of LinBox namespace
#endif // __DET_H
