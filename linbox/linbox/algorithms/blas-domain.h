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
#include <linbox/blackbox/permutation.h>
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
		Element _MOne;

	public:

		// Constructor of LSP.
		// Initialize the matrix L,S and the permutation P and some constant
		BlasMatrixDomain (const Field& F ) : _F(F) { F.init(_One,1UL); F.init(_Zero,0UL);F.init(_MOne,-1L);}
	    
		// Copy constructor
		BlasMatrixDomain (const BlasMatrixDomain<Field> & BMD) : _F(BMD._F), _One(BMD._One), _Zero(BMD._Zero), _MOne(BMD._MOne) {}


		// Field accessor
		Field& field() {return _F;}

			
		/*
		 * Basics operation available for BlasMatrix
		 */
 
		// multiplication
		// C = A*B
		template <class Operand, class Matrix>
		Operand& mul(Operand& C, const Matrix& A, const Operand& B) const;


		// multiplication with scaling
		// C = alpha.A*B
		template <class Operand, class Matrix>
		Operand& mul(Operand& C, const Element& alpha, const Matrix& A, const Operand& B) const;

		// axpy
		// D = C + A*B
		template <class Operand, class Matrix>
		Operand& axpy(Operand& D, const Matrix& A, const Operand& B, const Operand& C) const;

		// axpyin
		// C += A*B
		template <class Operand, class Matrix>
		Operand& axpyin(Operand& C, const Matrix& A, const Operand& B) const;
 
		// axmy
		// D= C - A*B
		template <class Operand, class Matrix>
		Operand& axmy(Operand& D, const Matrix& A, const Operand& B, const Operand& C) const;

		// axmyin
		// C-= A*B
		template <class Operand, class Matrix>
		Operand& axmyin(Operand& C, const Matrix& A, const Operand& B) const;
		
		//  general matrix-matrix multiplication and addition with scaling
		// D= beta.C + alpha.A*B
		template <class Operand, class Matrix>
		Operand& muladd(Operand& D, 
					   const Element& beta, const Operand& C,
					   const Element& alpha, const Matrix& A, const Operand& B) const;
		
		// C= beta.C + alpha.A*B
		template <class Operand, class Matrix>
		Operand& muladdin(const Element& beta, Operand& C,
					     const Element& alpha, const Matrix& A, const Operand& B) const;



		/*
		 * Solutions available for BlasMatrix 
		 */	

		// Inversion
		template <class Matrix>
		Matrix& inv( Matrix& Ainv, const Matrix& A) const;

		// Rank
		template <class Matrix>
		unsigned int rank(const Matrix& A) const;

		// in-place Rank (the matrix is modified)
		template <class Matrix>
		unsigned int rankin(Matrix& A) const;

		// determinant
		template <class Matrix>
		Element& det(const Matrix& A) const;

		//in-place Determinant (the matrix is modified)
		template <class Matrix>
		Element& detin(Matrix& A) const;
		
		/*
		 * Solvers for Matrix (respecting BlasMatrix interface) 
		 * with right or left Operand hand side
		 */ 

		// non-singular linear solve with matrix right hand side 
		// AX=B
		template <class Operand, class Matrix>
		Operand& left_solve (Operand& X, const Matrix& A, const Operand& B) const;
		
		// non-singular linear solve with matrix right hand side, the result is stored in-place in B
		// AX=B , (B<-X)
		template <class Operand,class Matrix>
		Operand& left_solve (const Matrix& A, Operand& B) const;
		
		// non-singular linear solve with matrix right hand side 
		// XA=B
		template <class Operand, class Matrix>
		Operand& right_solve (Operand& X, const Matrix& A, const Operand& B) const;
		
		// non-singular linear solve with matrix right hand side, the result is stored in-place in B
		// XA=B , (B<-X)
		template <class Operand, class Matrix>
		Operand& right_solve (const Matrix& A, Operand& B) const;
	       			
		
		/*
		 *  Method to apply Permutation
		 */
		// Apply a BlasPermutation matrix P to a dense matrix A: 
		// B = A.P 
		template <class Operand>
		Operand& applyRight(  Operand& B, const Operand& A, const BlasPermutation& P);

		// B = A.P^t
                template <class Operand>
		Operand& applyRightTranspose(  Operand& B, const Operand& A, const BlasPermutation& P);

		// B = P.A 
		template <class Operand>
		Operand& applyLeft(  Operand& B, const Operand& A, const BlasPermutation& P);
		
		// B = A.P^t
                template <class Operand>
		Operand& applyLeftTranspose(  Operand& B, const Operand& A, const BlasPermutation& P);
		
		// In place apply.
		// A = A.P 
		template <class Operand>
		Operand& applyinRight( Operand& A, const BlasPermutation& P);
		
		// A = A.P^t
                template <class Operand>
		Operand& applyinRightTranspose( Operand& A, const BlasPermutation& P);       

		// A = P.A 
		template <class Operand>
		Operand& applyinLeft( Operand& A, const BlasPermutation& P);
		
		// A = A.P^t
                template <class Operand>
		Operand& applyinLeftTranspose( Operand& A, const BlasPermutation& P);

		// Conversion from BlasPermutation to BlackBoxPermutation 
		//Permutation& convert ( Permutation& P, const BlasPermutation& BP );
		
	}; /* end of class BlasMatrixDomain */


} /* end of namespace LinBox */

#include <linbox/algorithms/blas-domain.inl>

#endif /* __BLAS_DOMAIN_H*/

