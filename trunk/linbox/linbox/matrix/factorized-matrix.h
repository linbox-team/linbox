/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/matrix/factorized-matrix.h
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


#ifndef __FACTORIZED_MATRIX_H
#define __FACTORIZED_MATRIX_H

#include <vector>

#include <linbox/matrix/blas-matrix.h>
#include <linbox/algorithms/blas-domain.h>
#include <linbox/fflapack/fflapack.h>

namespace LinBox{

	template <class Field,class Matrix>
	class LQUPMatrix {

	public:
		typedef typename Field::Element Element;
		typedef std::vector<size_t> BlasPermutation;

	protected:

		Field                   _F;
		BlasMatrix<Matrix>    &_LU;
		BlasPermutation         _P;
		BlasPermutation         _Q; 
		size_t                  _m;
		size_t                  _n;
		size_t               _rank;
    
	public:


		// Contruction of LQUP factorization of A (making a copy of A)
		LQUPMatrix (const Field& F, const BlasMatrix<Matrix>& A)
			: _F(F), _LU(*(new BlasMatrix<Matrix> (A))) , P(A.coldim()), Q(A.rowdim()), _m(A.rowdim()), _n(A.coldim())  {

			_rank= FFLAPACK::LUdivine( _F,FFLAS::FflasNonUnit, m, n, _LU.getPointer(),_LU.getStride(), &_P[0], FFLAPACK::FflapackLQUP, &_Q[0] );
			
		}

		// Contruction of LQUP factorization of A (in-place in A)
		LQUPMatrix (const Field& F, BlasMatrix<Matrix>& A)
			: _F(F), _LU(A) , P(A.coldim()), Q(A.rowdim()), _m(A.rowdim()), _n(A.coldim())  {

			_rank= FFLAPACK::LUdivine( _F,FFLAS::FflasNonUnit, m, n, _LU.getPointer(),_LU.getStride(), &_P[0], FFLAPACK::FflapackLQUP, &_Q[0] );
			
		}

		// get the field on which the factorization is done
		Field& field() {return _F;}    

		// get row dimension
		size_t rowdim() {return _m;}

		// get column dimension
		size_t coldim() const {return _n;}
    
		// get the rank of matrix
		size_t getrank() const {return _rank;}
    
		// get the permutation P
		const BlasPermutation& getP() const {return _P;}
       
		// get the permutation Q
		const BlasPermutation& getQ() const  {return _Q;}

		// get the Matrix L
		const BlasMatrix<Matrix>& getL() const;

		// get the matrix U
		const BlasMatrix<Matrix>& getU() const;

		// get the matrix S (from the LSP factorization of A deduced from LQUP)
		const BlasMatrix<Matrix>& getS() const;



		/*
		 * Solvers with Matrices
		 */
		// solve AX=B
		bool left_solve(BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& B) const;

		// solve AX=B (X is stored in B)
		bool left_solve(BlasMatrix<Matrix>& B) const;

		// solve XA=B
		bool right_solve(BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& B) const;

		// solve XA=B (X is stored in B)
		bool right_solve(BlasMatrix<Matrix>& B) const;


		// solve LX=B (L from LQUP)
		bool left_Lsolve(BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& B) const;
		
		// solve LX=B (L from LQUP) (X is stored in B)
		bool left_Lsolve(BlasMatrix<Matrix>& B) const;

		// solve XL=B (L from LQUP)
		bool right_Lsolve(BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& B) const;
		
		// solve XL=B (L from LQUP) (X is stored in B)
		bool right_Lsolve(BlasMatrix<Matrix>& B) const;


		// solve UX=B (U from LQUP)
		bool left_Usolve(BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& B) const;
		
		// solve UX=B (U from LQUP) (X is stored in B)
		bool rleft_Usolve(BlasMatrix<Matrix>& B) const;

		// solve XU=B (U from LQUP)
		bool right_Usolve(BlasMatrix<Matrix>& X, const BlasMatrix<Matrix>& B) const;
		
		// solve XU=B (U from LQUP) (X is stored in B)
		bool right_Usolve(BlasMatrix<Matrix>& B) const;


		/*
		 * Solvers with vectors
		 */
		
		// solve Ax=b
		bool left_solve(std::vector<Element>& x, const std::vector<Element>& b) const;
		
		// solve Ax=b (x is stored in b)
		bool left_solve(std::vector<Element>& b) const;

		// solve xA=b
		bool right_solve(std::vector<Element>& x, const std::vector<Element>& b) const;

		// solve xA=b (x is stored in b)
		bool right_solve(std::vector<Element>& b) const;


		// solve Lx=b (L from LQUP)
		bool left_Lsolve(std::vector<Element>& x, const std::vector<Element>& b) const;
		
		// solve Lx=b (L from LQUP) (x is stored in b)
		bool left_Lsolve(std::vector<Element>& b) const;

		// solve xL=b (L from LQUP)
		bool right_Lsolve(std::vector<Element>& x, const std::vector<Element>& b) const;
		
		// solve xL=b (L from LQUP) (x is stored in b)
		bool right_Lsolve(std::vector<Element>& b) const;


		// solve Ux=b (U from LQUP)
		bool left_Usolve(std::vector<Element>& x, const std::vector<Element>& b) const;
		
		// solve Ux=b (U from LQUP) (x is stored in b)
		bool rleft_Usolve(std::vector<Element>& b) const;

		// solve xU=b (U from LQUP)
		bool right_Usolve(std::vector<Element>& x, const std::vector<Element>& b) const;
		
		// solve xU=b (U from LQUP) (x is stored in b)
		bool right_Usolve(std::vector<Element>& b) const;


	}; //

} //end of namespace LinBox

#include <linbox/matrix/matrix-factorized.inl>

#endif
