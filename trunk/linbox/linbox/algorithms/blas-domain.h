/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/blas-domain.h
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@ens-lyon.fr
 *               Clément Pernet clement.pernet@imag.fr
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

#ifndef __BLAS_MATRIX_DOMAIN_H
#define __BLAS_MATRIX_DOMAIN_H

#include <iostream>
#include <vector>
#include <linbox/fflapack/fflapack.h>
#include <linbox/fflas/fflas.h>
#include <linbox/blackbox/permutation;h>
#include <linbox/matrix/blas-matrix.h>






namespace LinBox {

	template <class Field>
	class BlasMatrixDomain {

	public:
		typedef typename Field::Element         Element;
		typedef std::vector<size_t>     BlasPermutation;

	private:
    
		Field  _F;

	public:

		// Constructor of LSP.
		// Initialize the matrix L,S and the permutation P and some constant
		BlasMatrixDomain (const Field& F ) : _F(F) {}
	    
		// Copy constructor
		BlasMatrixDomain (const BlasMatrixDomain<Field> & BMD) : _F(BMD._F) {}


		// Field accessor
		Field& field() {return _F;}

		// LSP Factorization 
		
		const unsigned int LSP (const BlasMatrix<Matrix>& A, BlasMatrix<Matrix>& L, BlasMatrix<Matrix>& S, BlasPermutation& P) const;					
		// LQUP Factorization
		const unsigned int LQUP (const BlasMatrix<Matrix>& A, BlasMatrix<Matrix>& L, BlasPermutation& Q, BlasMatrix<Matrix>& U, BlasPermutation& P) const;
		
		// in-place LQUP Factorization (L is in compressed format)
		const unsigned int LQUPin (BlasMatrix<Matrix>& A, BlasPermutation& Q, BlasPermutation& P) const;

		// Inversion
		const BlasMatrix<Matrix>& inv(const BlasMatrix<Matrix>& A, BlasMatrix<Matrix>& Ainv) const;

		// Rank
		const unsigned int rank(const BlasMatrix<Matrix>& A) const;

		// in-place Rank (the matrix is modified)
		const unsigned int rankin(BlasMatrix<Matrix>& A) const;

		// determinant
		const Element& det(const BlasMatrix<Matrix>& A) const;

		//in-place Determinant (the matrix is modified)
		const Element& detin(BlasMatrix<Matrix>& A) const;
		
		// non-singular linear solve with matrix right hand side 
		bool solve (BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& A, const BlasMatrix<Matrix>& B) const;

		// non-singular linear solve with vector right hand side
		bool solve (std::vector<Matrix>& x, const BlasMatrix<Matrix>& A, const std::vector<Matrix& b) const;

		// non-singular linear solve with matrix right hand side, the result is stored in-place in B
		bool solve (const BlasMatrix<Matrix>& A, const BlasMatrix<Matrix>& B) const;

		// non-singular linear solve with vector right hand side
		bool solve (const BlasMatrix<Matrix>& A, const std::vector<Matrix& b) const;

		// Apply a BlasPermutation matrix P to a dense matrix A: 
		// B = A.P 
		BlasMatrix<Matrix>& applyRight(  BlasMatrix<Matrix>& B, const BlasMatrix<Matrix>& A, const BlasPermutation& P);

		// B = A.P^t
                BlasMatrix<Matrix>& applyRightTranspose(  BlasMatrix<Matrix>& B, const BlasMatrix<Matrix>& A, const BlasPermutation& P);

		// B = P.A 
		BlasMatrix<Matrix>& applyLeft(  BlasMatrix<Matrix>& B, const BlasMatrix<Matrix>& A, const BlasPermutation& P);
		
		// B = A.P^t
                BlasMatrix<Matrix>& applyLeftTranspose(  BlasMatrix<Matrix>& B, const BlasMatrix<Matrix>& A, const BlasPermutation& P);
		
		// In place apply.
		// A = A.P 
		BlasMatrix<Matrix>& applyinRight( BlasMatrix<Matrix>& A, const BlasPermutation& P);
		
		// A = A.P^t
                BlasMatrix<Matrix>& applyinRightTranspose( BlasMatrix<Matrix>& A, const BlasPermutation& P);       

		// A = P.A 
		BlasMatrix<Matrix>& applyinLeft( BlasMatrix<Matrix>& A, const BlasPermutation& P);
		
		// A = A.P^t
                BlasMatrix<Matrix>& applyinLeftTranspose( BlasMatrix<Matrix>& A, const BlasPermutation& P);

		// Conversion from BlasPermutation to BlackBoxPermutation 
		Permutation& convert ( Permutation& P, const BlasPermutation& BP );
		
	}; /* end of class BlasMatrixDomain */


} /* end of namespace LinBox */

#endif /* __BLAS_DOMAIN_H*/
