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
#include <linbox/util/debug.h>





namespace LinBox {
		

	template <class Field>
	class BlasMatrixDomain {

	public:
		typedef typename Field::Element         Element;
		typedef std::vector<size_t>     BlasPermutation;

	private:
    
		Field  _F;
		Element _One;
		Element _Zero;		

	public:

		// Constructor of LSP.
		// Initialize the matrix L,S and the permutation P and some constant
		BlasMatrixDomain (const Field& F ) : _F(F) { F.init(_One,1UL); F.init(_Zero,0UL);}
	    
		// Copy constructor
		BlasMatrixDomain (const BlasMatrixDomain<Field> & BMD) : _F(BMD._F), _One(BMD._One), _Zero(BMD._Zero) {}


		// Field accessor
		Field& field() {return _F;}

			
		/*
		 * Solutions available for BlasMatrix 
		 */	

		// Inversion
		template <class Matrix>
		const BlasMatrix<Matrix>& inv(const BlasMatrix<Matrix>& A, BlasMatrix<Matrix>& Ainv) const;

		// Rank
		template <class Matrix>
		const unsigned int rank(const BlasMatrix<Matrix>& A) const;

		// in-place Rank (the matrix is modified)
		template <class Matrix>
		const unsigned int rankin(BlasMatrix<Matrix>& A) const;

		// determinant
		template <class Matrix>
		const Element& det(const BlasMatrix<Matrix>& A) const;

		//in-place Determinant (the matrix is modified)
		const Element& detin(BlasMatrix<Matrix>& A) const;
		
		/*
		 * Solvers with Matrix right or left hand side
		 */ 
		// non-singular linear solve with matrix right hand side 
		// AX=B
		template <class Matrix>
		bool left_solve (BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& A, const BlasMatrix<Matrix>& B) const;
		
		// non-singular linear solve with matrix right hand side, the result is stored in-place in B
		// AX=B , (B<-X)
		template <class Matrix>
		bool left_solve (const BlasMatrix<Matrix>& A, const BlasMatrix<Matrix>& B) const;
		
		// non-singular linear solve with matrix right hand side 
		// XA=B
		template <class Matrix>
		bool right_solve (BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& A, const BlasMatrix<Matrix>& B) const;
		
		// non-singular linear solve with matrix right hand side, the result is stored in-place in B
		// XA=B , (B<-X)
		template <class Matrix>
		bool right_solve (const BlasMatrix<Matrix>& A, const BlasMatrix<Matrix>& B) const;

		/*
		 * Solvers with vectors right or left hand side
		 */
		// non-singular linear solve with matrix right hand side 
		// Ax=b
		template <class Matrix>
		bool left_solve (std::vector<Element>& X, const BlasMatrix<Matrix>& A, const std::vector<Element>& B) const;
		
		// non-singular linear solve with matrix right hand side, the result is stored in-place in B
		// Ax=b , (b<-x)
		template <class Matrix>
		bool left_solve (const BlasMatrix<Matrix>& A, const std::vector<Element>& B) const;
		
		// non-singular linear solve with matrix right hand side 
		// xA=b
		template <class Matrix>
		bool right_solve (std::vector<Element>& X, const BlasMatrix<Matrix>& A, const std::vector<Element>& B) const;
		
		// non-singular linear solve with matrix right hand side, the result is stored in-place in B
		// xA=b , (b<-x)
		template <class Matrix>
		bool right_solve (const BlasMatrix<Matrix>& A, const std::vector<Element>& B) const;


		
		// Solvers available for Triangular Blas Matrix 

		/*
		 * with Matrix right or left hand side
		 */ 
		// non-singular linear solve with matrix right hand side 
		// TX=B
		template <class Matrix>
		bool left_solve (BlasMatrix<Matrix>& X, const TriangularBlasMatrix<Matrix>& T, const BlasMatrix<Matrix>& B) const;
		
		// non-singular linear solve with matrix right hand side, the result is stored in-place in B
		// TX=B , (B<-X)
		template <class Matrix>
		bool left_solve (const TriangularBlasMatrix<Matrix>& T, const BlasMatrix<Matrix>& B) const;
		
		// non-singular linear solve with matrix right hand side 
		// XT=B
		template <class Matrix>
		bool right_solve (BlasMatrix<Matrix>& X, const TriangularBlasMatrix<Matrix>& T, const BlasMatrix<Matrix>& B) const;
		
		// non-singular linear solve with matrix right hand side, the result is stored in-place in B
		// XT=B , (B<-X)
		template <class Matrix>
		bool right_solve (const TriangularBlasMatrix<Matrix>& T, const BlasMatrix<Matrix>& B) const;

		/*
		 * with vectors right or left hand side
		 */
		// non-singular linear solve with matrix right hand side 
		// Tx=b
		template <class Matrix>
		bool left_solve (std::vector<Element>& x, const TriangularBlasMatrix<Matrix>& T, const std::vector<Element>& b) const;
		
		// non-singular linear solve with matrix right hand side, the result is stored in-place in b
		// Tx=b , (b<-x)
		template <class Matrix>
		bool left_solve (const TriangularBlasMatrix<Matrix>& T, const std::vector<Element>& b) const;
		
		// non-singular linear solve with matrix right hand side 
		// xT=b
		template <class Matrix>
		bool right_solve (std::vector<Element>& x, const TriangularBlasMatrix<Matrix>& T, const std::vector<Element>& b) const;
		
		// non-singular linear solve with matrix right hand side, the result is stored in-place in b
		// xT=b , (b<-x)
		template <class Matrix>
		bool right_solve (const TriangularBlasMatrix<Matrix>& T, const std::vector<Element>& b) const;

	       
		


		/*
		 *  Method to apply Permutation
		 */
		// Apply a BlasPermutation matrix P to a dense matrix A: 
		// B = A.P 
		template <class Matrix>
		BlasMatrix<Matrix>& applyRight(  BlasMatrix<Matrix>& B, const BlasMatrix<Matrix>& A, const BlasPermutation& P);

		// B = A.P^t
                template <class Matrix>
		BlasMatrix<Matrix>& applyRightTranspose(  BlasMatrix<Matrix>& B, const BlasMatrix<Matrix>& A, const BlasPermutation& P);

		// B = P.A 
		template <class Matrix>
		BlasMatrix<Matrix>& applyLeft(  BlasMatrix<Matrix>& B, const BlasMatrix<Matrix>& A, const BlasPermutation& P);
		
		// B = A.P^t
                template <class Matrix>
		BlasMatrix<Matrix>& applyLeftTranspose(  BlasMatrix<Matrix>& B, const BlasMatrix<Matrix>& A, const BlasPermutation& P);
		
		// In place apply.
		// A = A.P 
		template <class Matrix>
		BlasMatrix<Matrix>& applyinRight( BlasMatrix<Matrix>& A, const BlasPermutation& P);
		
		// A = A.P^t
                template <class Matrix>
		BlasMatrix<Matrix>& applyinRightTranspose( BlasMatrix<Matrix>& A, const BlasPermutation& P);       

		// A = P.A 
		template <class Matrix>
		BlasMatrix<Matrix>& applyinLeft( BlasMatrix<Matrix>& A, const BlasPermutation& P);
		
		// A = A.P^t
                template <class Matrix>
		BlasMatrix<Matrix>& applyinLeftTranspose( BlasMatrix<Matrix>& A, const BlasPermutation& P);

		// Conversion from BlasPermutation to BlackBoxPermutation 
		Permutation& convert ( Permutation& P, const BlasPermutation& BP );
		
	}; /* end of class BlasMatrixDomain */


} /* end of namespace LinBox */

#endif /* __BLAS_DOMAIN_H*/
