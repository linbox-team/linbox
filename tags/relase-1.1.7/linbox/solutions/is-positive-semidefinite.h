/* Copyright (C) 2010 LinBox
 *
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

#ifndef __LINBOX_is_positive_semidefinite_H
#define __LINBOX_is_positive_semidefinite_H

#include "linbox/util/error.h"
#include <linbox/algorithms/matrix-hom.h>
#include "linbox/algorithms/signature.h"

namespace LinBox
{
	// for specialization with respect to the DomainCategory
    template< class Blackbox, class isPositiveSemiDefiniteMethod, class DomainCategory>
    bool isPositiveSemiDefinite (
		const Blackbox        &A,
		const DomainCategory  &tag,
		const isPositiveSemiDefiniteMethod  &M);

	/** Compute the isPositiveSemiDefinite of A
	 *
	 * The isPositiveSemiDefinite of a linear operator A, represented as a
	 * black box, is computed over the ring or field of A.
	 *
	 * @param r OUTPUT instance into which to store the result r
	 * @param A Black box of which to compute the isPositiveSemiDefinite
	 * @param M may be a Method::Hybrid (SemiDefault), Method::Blackbox, Method::Elimination, or of other method type.
         \ingroup isPositiveSemiDefinites
        */
    template <class Blackbox, class MyMethod>
    bool isPositiveSemiDefinite (
		const Blackbox                              &A,
		const MyMethod                           &M) 
    {
        return isPositiveSemiDefinite( A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
    }

	// The isPositiveSemiDefinite with SemiDefault Method 
    template<class Blackbox>
    bool isPositiveSemiDefinite ( const Blackbox  &A) {
        return isPositiveSemiDefinite(A, 
		Method::Hybrid());
    }

	// The isPositiveSemiDefinite for ModularTag (is nonsense)
    template<class Blackbox, class MyMethod>
    bool isPositiveSemiDefinite (
        const Blackbox                            &A,
        const RingCategories::ModularTag          &tag,
		const MyMethod& M)
    {
		//commentator << "nonsense!!"
		throw (LinboxError("isPositiveSemiDefinite: Integer matrix required"));
        return false;
    }

	// The isPositiveSemiDefinite with Hybrid Method 
    template<class Blackbox>
    bool isPositiveSemiDefinite (
        const Blackbox 			&A,
        const RingCategories::IntegerTag          &tag,
		const Method::Hybrid& M)
    {
		// should try a modular minpoly and decide on the degree of that...
        if (A.rowdim() != A.coldim()) return false;
		// this crude size check can be refined
		if (A.coldim() > 7000) return isPositiveSemiDefinite(A, tag, Method::Blackbox(M));
		else return isPositiveSemiDefinite(A, tag, Method::Elimination(M));
    }

	// The isPositiveSemiDefinite with Elimination Method 
    template<class Blackbox>
    bool isPositiveSemiDefinite (
		const Blackbox                            &A,
		const RingCategories::IntegerTag          &tag,
		const Method::Elimination& M)
    {
		// this can be a hybrid of EliminationMinpoly and BlasElimination (which means use LU here)
		// It will be faster to do EliminationMinpoly when deg(m_A) is low.

		// right now it is just BlasElimination
        return isPositiveSemiDefinite(A, tag, Method::BlasElimination(M));
    }

	// The isPositiveSemiDefinite with BlackBox Method 
    template<class Blackbox>
    bool isPositiveSemiDefinite (
		const Blackbox                      &A,
		const RingCategories::IntegerTag    &tag,
		const Method::Blackbox              &M)
    {
        return isPositiveSemiDefinite(A, tag, Method::Wiedemann(M));
    }


	// The isPositiveSemiDefinite with Wiedemann, finite field.
    template <class Blackbox>
    bool isPositiveSemiDefinite (
		const Blackbox                      &A,
		const RingCategories::IntegerTag    &tag,
		const Method::Wiedemann             &M)
    {
		// call Wiedemann code
		return Signature::isPosSemiDef(A, Signature::Minpoly_Method() );
	}

	// the isPositiveSemiDefinite with Blas. 
    template <class Blackbox>
    bool isPositiveSemiDefinite (
		const Blackbox                      &A,
		const RingCategories::IntegerTag    &tag,
		const Method::BlasElimination       &M)
    {
		// call BlasElimination code
		DenseMatrix<typename Blackbox::Field> DA(A.field(), A.rowdim(), A.coldim());
		MatrixHom::map(DA, A, A. field());
		return Signature::isPosSemiDef(DA, Signature::BLAS_LPM_Method() );
	}
	
	// the isPositiveSemiDefinite with Blas, DenseMatrix
    template <class Ring>
    bool isPositiveSemiDefinite (
		const DenseMatrix<Ring> &A,
		const RingCategories::IntegerTag    &tag,
		const Method::BlasElimination       &M)
    {
		// call BlasElimination code
		return Signature::isPosSemiDef(A, Signature::BLAS_LPM_Method() );
	}
	
} // end of LinBox namespace
#endif // __LINBOX_is_positive_semidefinite_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
