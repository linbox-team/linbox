/* Copyright (C) 2010 LinBox
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

#ifndef __LINBOX_is_positive_definite_H
#define __LINBOX_is_positive_definite_H

#include "linbox/util/error.h"
#include <linbox/algorithms/matrix-hom.h>
#include "linbox/algorithms/signature.h"

namespace LinBox
{
	// for specialization with respect to the DomainCategory
    template< class Blackbox, class isPositiveDefiniteMethod, class DomainCategory>
    bool isPositiveDefinite (
		const Blackbox        &A,
		const DomainCategory  &tag,
		const isPositiveDefiniteMethod  &M);

	/** Compute the isPositiveDefinite of A
	 *
	 * The isPositiveDefinite of a linear operator A, represented as a
	 * black box, is computed over the ring or field of A.
	 *
	 * @param r OUTPUT instance into which to store the result r
	 * @param A Black box of which to compute the isPositiveDefinite
	 * @param M may be a Method::Hybrid (default), Method::Blackbox, Method::Elimination, or of other method type.
         \ingroup isPositiveDefinites
        */
    template <class Blackbox, class MyMethod>
    bool isPositiveDefinite (
		const Blackbox                              &A,
		const MyMethod                           &M) 
    {
        return isPositiveDefinite( A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
    }

	// The isPositiveDefinite with default Method 
    template<class Blackbox>
    bool isPositiveDefinite ( const Blackbox  &A) {
        return isPositiveDefinite(A, 
		Method::Hybrid());
    }

	// The isPositiveDefinite for ModularTag (is nonsense)
    template<class Blackbox, class MyMethod>
    bool isPositiveDefinite (
        const Blackbox                            &A,
        const RingCategories::ModularTag          &tag,
		const MyMethod& M)
    {
		//commentator << "nonsense!!"
		throw (LinboxError("isPositiveDefinite: Integer matrix required"));
        return false;
    }

	// The isPositiveDefinite with Hybrid Method 
    template<class Blackbox>
    bool isPositiveDefinite (
        const Blackbox 			&A,
        const RingCategories::IntegerTag          &tag,
		const Method::Hybrid& M)
    {
		// should try a modular minpoly and decide on the degree of that...
        if (A.rowdim() != A.coldim()) return false;
		// this crude size check can be refined
		if (A.coldim() > 7000) return isPositiveDefinite(A, tag, Method::Blackbox(M));
		else return isPositiveDefinite(A, tag, Method::Elimination(M));
    }

	// The isPositiveDefinite with Elimination Method 
    template<class Blackbox>
    bool isPositiveDefinite (
		const Blackbox                            &A,
		const RingCategories::IntegerTag          &tag,
		const Method::Elimination& M)
    {
		// this can be a hybrid of EliminationMinpoly and BlasElimination (which means use LU here)
		// It will be faster to do EliminationMinpoly when deg(m_A) is low.

		// right now it is just BlasElimination
        return isPositiveDefinite(A, tag, Method::BlasElimination(M));
    }

	// The isPositiveDefinite with BlackBox Method 
    template<class Blackbox>
    bool isPositiveDefinite (
		const Blackbox                      &A,
		const RingCategories::IntegerTag    &tag,
		const Method::Blackbox              &M)
    {
        return isPositiveDefinite(A, tag, Method::Wiedemann(M));
    }


	// The isPositiveDefinite with Wiedemann, finite field.
    template <class Blackbox>
    bool isPositiveDefinite (
		const Blackbox                      &A,
		const RingCategories::IntegerTag    &tag,
		const Method::Wiedemann             &M)
    {
		// call Wiedemann code
		return Signature::isPosDef(A, Signature::Minpoly_Method() );
	}

	// the isPositiveDefinite with Blas. 
    template <class Blackbox>
    bool isPositiveDefinite (
		const Blackbox                      &A,
		const RingCategories::IntegerTag    &tag,
		const Method::BlasElimination       &M)
    {
		// call BlasElimination code
		DenseMatrix<typename Blackbox::Field> DA(A.field(), A.rowdim(), A.coldim());
		MatrixHom::map(DA, A, A. field());
		return Signature::isPosDef(DA, Signature::BLAS_LPM_Method() );
    }
	
	// the isPositiveDefinite with Blas, DenseMatrix
    template <class Ring>
    bool isPositiveDefinite (
		const DenseMatrix<Ring> &A,
		const RingCategories::IntegerTag    &tag,
		const Method::BlasElimination       &M)
    {
		// call BlasElimination code
		return Signature::isPosDef(A, Signature::BLAS_LPM_Method() );
	}
	
} // end of LinBox namespace

#endif // __LINBOX_is_positive_definite_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
