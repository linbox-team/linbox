/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/solutions/minpoly.h
 * Copyright (C) 1999, 2001 Jean-Guillaume Dumas, Bradford Hovinen
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __MINPOLY_H
#define __MINPOLY_H

#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/blackbox-container-symmetric.h"
#include "linbox/algorithms/massey-domain.h"     // massey recurring sequence solver
#include "linbox/solutions/methods.h"
#include "linbox/util/commentator.h"
#include <linbox/algorithms/minpoly-integer.h>


namespace LinBox 
{
	
	/** @memo Minimal polynomial of a blackbox linear operator A.
	 * @doc The resulting polynomial is a vector of coefficients.
	 * Somewhere we should document our handling of polys.
	 */
	template < class Blackbox, class Polynomial, class FieldCategoryTag>
	Polynomial &minpoly (Polynomial& P,
			     const Blackbox& A,
			     const FieldCategoryTag& tag,
			     const Method::Wiedemann& M = Method::Wiedemann ());
	
	template < class Blackbox, class Polynomial>
	Polynomial &minpoly (Polynomial& P,
			     const Blackbox& A,
			     const Method::Wiedemann& M = Method::Wiedemann ()) {

		return minpoly (P, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
	}

	template < class Blackbox, class Polynomial>
	Polynomial &minpoly (Polynomial& P,
			     const Blackbox& A,
			     RingCategories::IntegerTag tag,
			     const Method::Wiedemann& M = Method::Wiedemann ())
	{	
		typedef Modular<double> ModularField;
		MinPoly<typename Blackbox::Field::Element, ModularField>::minPoly(P, A);

		return P;
	}

	template < class Blackbox, class Polynomial, class FieldCategoryTag>
	Polynomial &minpolySymmetric (Polynomial& P,
				      const Blackbox& A,
				      FieldCategoryTag tag,
				      const Method::Wiedemann& M = Method::Wiedemann ());

	template < class Blackbox, class Polynomial>
	Polynomial &minpolySymmetric (Polynomial& P,
				      const Blackbox& A,
				      const Method::Wiedemann& M = Method::Wiedemann ()) 
	{

		minpolySymmetric(P, A,  typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
		return P;
	}
				 
	template < class Blackbox, class Polynomial>
	Polynomial &minpolySymmetric (Polynomial& P,
				      const Blackbox& A,
				      RingCategories::IntegerTag tag,
				      const Method::Wiedemann& M = Method::Wiedemann ())
	{	
		typedef typename Blackbox::Field::Element Integer;
		typedef Modular<double> ModularField;
		MinPoly<Integer, ModularField>::minPoly(P, A);

		return P;
	}

	template < class Blackbox, class Polynomial>
	Polynomial &minpoly (Polynomial& P,
			     const Blackbox& A,
			     RingCategories::ModularTag tag,
			     const Method::Wiedemann& M = Method::Wiedemann ())
	{
		typedef typename Blackbox::Field Field;
		typename Field::RandIter i (A.field());
		unsigned long            deg;

		commentator.start ("Minimal polynomial", "minpoly");

		BlackboxContainer<Field, Blackbox> TF (&A, A.field(), i);
		MasseyDomain< Field, BlackboxContainer<Field, Blackbox> > WD (&TF, M.earlyTermThreshold ());

		WD.minpoly (P, deg);

#ifdef INCLUDE_TIMING
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for applies:      " << TF.applyTime () << endl;
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for dot products: " << TF.dotTime () << endl;
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for discrepency:  " << WD.discrepencyTime () << endl;
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for LSR fix:      " << WD.fixTime () << endl;
#endif // INCLUDE_TIMING

		commentator.stop ("done", NULL, "minpoly");

		return P;
	}

	template < class Blackbox, class Polynomial>
	Polynomial &minpolySymmetric (Polynomial& P,
				      const Blackbox& A,
				      RingCategories::ModularTag tag,
				      const Method::Wiedemann& M = Method::Wiedemann ())
	{
		typedef typename Blackbox::Field Field;
		typename Field::RandIter i (A.field());
		unsigned long            deg;

		commentator.start ("Minimal polynomial", "minpoly");

		BlackboxContainerSymmetric<Field, Blackbox> TF (&A, A.field(), i);
		MasseyDomain< Field, BlackboxContainerSymmetric<Field, Blackbox> > WD (&TF, M.earlyTermThreshold ());

		WD.minpoly (P, deg);

#ifdef INCLUDE_TIMING
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for applies:      " << TF.applyTime () << endl;
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for dot products: " << TF.dotTime () << endl;
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for discrepency:  " << WD.discrepencyTime () << endl;
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for LSR fix:      " << WD.fixTime () << endl;
#endif // INCLUDE_TIMING

		commentator.stop ("done", NULL, "minpoly");

		return P;
	}

}
#endif // __MINPOLY_H
