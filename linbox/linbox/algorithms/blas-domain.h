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

/** This class provide the decomposition M = LSP of a dense matrix which is stored contiguously by row
 * and where row size is lower or equal than the column size. 
 *
 * The class is templatized by a Matrix type which must have the function FullIterator(). This function returns a pointer
 * to the 1st element of the matrix.
 */

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
		BlasMatrixDomain (const BlasMatrixDomain & BMD) : _F(BMD._F) {}


		// Field accessor
		Field& field() {return _F;}

		// LSP Factorization 
		const unsigned int LSP (const BlasMatrix<Element>& A, BlasMatrix<Element>& L, BlasMatrix<Element>& S, BlasPermutation& P) const;					
		// in-place LSP Factorization (L is in compressed format)
		const unsigned int LSPin (BlasMatrix<Element>& A, BlasPermutation& P) const;

		// LQUP Factorization
		const unsigned int LQUP (const BlasMatrix<Element>& A, BlasMatrix<Element>& L, BlasPermutation& Q, BlasMatrix<Element>& U, BlasPermutation& P) const;
		
		// in-place LQUP Factorization (L is in compressed format)
		const unsigned int LQUPin (BlasMatrix<Element>& A, BlasPermutation& Q, BlasPermutation& P) const;

		// Inversion
		const BlasMatrix<Element>& Inversion(const BlasMatrix<Element>& A, BlasMatrix<Element>& Ainv) const;


		// Rank
		const unsigned int rank(const BlasMatrix<Element>& A) const;
		
		// in-place Rank (the matrix is modified)
		const unsigned int rankin(BlasMatrix<Element>& A) const;

		// determinant
		const Element& det(const BlasMatrix<Element>& A) const;

		//in-place determinant (the matrix is modified)
		const Element& detin(BlasMatrix<Element>& A) const;


		


	}; /* end of class BlasMatrixDomain */


} /* end of namespace LinBox */

#endif /* __BLAS_DOMAIN_H*/
