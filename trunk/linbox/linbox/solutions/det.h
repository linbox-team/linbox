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
#include "linbox/solutions/methods.h"

#include "linbox/blackbox/blas-blackbox.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h"
#include "linbox/algorithms/cra.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/prime-stream.h"
#include "linbox/util/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
	// Methods for rank
	struct RankMethod
	{
		typedef WiedemannTraits Wiedemann; // should allow for return of a probability figure.
		struct BlasElimination{};
	};
	// for specialization with respect to the DomainCategory
	template< class Blackbox, class RankMethod, class DomainCategory>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &d, 
						const Blackbox                              &A,
						const DomainCategory                      &tag,
						const RankMethod                          &M);

	/** Compute the determinant of A
	 *
	 * The determinant of a linear operator A, represented as a
	 * black box, is computed over the ring or field of A.
	 *
	 * @param d Field element into which to store the result
	 * @param A Black box of which to compute the determinant
	 * @param M may be a RankMethod::BlasElimination (default) or a RankMethod::Wiedemann.
	\ingroup solutions
	 */
	template <class Blackbox, class Method>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &d, 
						const Blackbox                              &A,				
						const Method                           &M) 
	{
		return det(d, A, FieldTraits<typename Blackbox::Field>::categoryTag(), M);
	}

	// The det with default Method 
	template<class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &d, 
						const Blackbox                               &A)
	{
		return det(d, A, FieldTraits<typename Blackbox::Field>(), RankMethod::BlasElimination());
	}

	// The det with Wiedemann, finite field.
	template <class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &d, 
						const Blackbox                              &A,
						const RingCategories::ModularTag          &tag,
						const RankMethod::Wiedemann                &M) 
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

		std::vector<typename Field::Element> diag (A.coldim ());

		typename Field::Element pi;
		size_t i;
		size_t iternum = 1;
		do {
			F.init (pi, 1);
			for (i = 0; i < A.coldim (); i++) {
				do iter.random (diag[i]); while (F.isZero (diag[i]));
				F.mulin (pi, diag[i]);
			}
	
			Diagonal<Field> D (F, diag);

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

		F.div (d, phi[0], pi);

		commentator.stop ("done", NULL, "det");

		return d;
	}



	// the det with Blas, finite field.
	template <class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &d,
						const Blackbox                              &A,
						const RingCategories::ModularTag          &tag,
						const RankMethod::BlasElimination           &M) 
	{
		typedef typename Blackbox::Field Field;
		Field F = A.field();
		
		commentator.start ("Determinant", "det");

		linbox_check (A.coldim () == A.rowdim ());

		BlasMatrix<typename Field::Element> B(A);
		BlasMatrixDomain<Field> BMD(F);
		d= BMD.det(B);
		commentator.stop ("done", NULL, "det");

		return d;
	}

	
	/// A will be modified.
	/** 
	 \todo This should work for a DenseMatrix.
	 \returns d determinant of A.
	 \param A this BlasBlackbox matrix will be modified in place in the process.
	\ingroup solutions
 	*/
	template <class Field>
	typename Field::Element &detin (typename Field::Element             &d,
					BlasBlackbox<Field>                   &A)
	{
		Field F = A.field();
		
		commentator.start ("Determinant", "det");
		linbox_check (A.coldim () == A.rowdim ());

		BlasMatrixDomain<Field> BMD(F);
		d= BMD.detin(static_cast<BlasMatrix<typename Field::Element>& > (A));
		commentator.stop ("done", NULL, "det");

		return d;
	}
} // end of LinBox namespace 

#include "linbox/algorithms/cra-det-integer.h"

namespace LinBox {
	
	template <class Blackbox, class Method>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &d,
						const Blackbox                              &A,
						const RingCategories::IntegerTag          &tag,
						const Method                           &M) 
	{
		return cra_det (d, A, M);
	}


} // end of LinBox namespace
#endif // __DET_H
