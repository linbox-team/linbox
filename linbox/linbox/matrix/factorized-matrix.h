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

	template <class Field, class Operand> 
	class FactorizedMatrixLeftSolve {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::element>& _LU, 
				  Operand& B ) const;
	}; // end of class FactorizedMatrixLeftSolve

	template <class Field, class Operand> 
	class FactorizedMatrixRightSolve {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::element>& _LU, 
				  Operand& B ) const;
	}; // end of class FactorizedMatrixRightSolve

	template <class Field, class Operand> 
	class FactorizedMatrixLeftLSolve {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::element>& _LU, 
				  Operand& B ) const;
	}; // end of class FactorizedMatrixLeftLSolve

	template <class Field, class Operand> 
	class FactorizedMatrixRightLSolve {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::element>& _LU, 
				  Operand& B ) const;
	}; // end of class FactorizedMatrixRightLsolve

	template <class Field, class Operand> 
	class FactorizedMatrixLeftUSolve {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::element>& _LU, 
				  Operand& B ) const;
	}; // end of class FactorizedMatrixLeftUSolve

	template <class Field, class Operand> 
	class FactorizedMatrixRightUSolve {
	public:
		bool operator() ( const Field& F, 
				  const BlasMatrix<typename Field::element>& _LU, 
				  Operand& B ) const;
	}; // end of class FactorizedMatrixRightUSolve
	
	template <class Field>
	class LQUPMatrix {

	public:
		typedef typename Field::Element Element;
		typedef std::vector<size_t> BlasPermutation;

	protected:

		Field                   _F;
		BlasMatrix<Element>    &_LU;
		BlasPermutation         _P;
		BlasPermutation         _Q; 
		size_t                  _m;
		size_t                  _n;
		size_t               _rank;
    
	public:


		// Contruction of LQUP factorization of A (making a copy of A)
		LQUPMatrix (const Field& F, const BlasMatrix<Element>& A)
			: _F(F), _LU(*(new BlasMatrix<Element> (A))) , 
			  P(A.coldim()), Q(A.rowdim()), _m(A.rowdim()), _n(A.coldim())  {

			_rank= FFLAPACK::LUdivine( _F,FFLAS::FflasNonUnit, m, n, 
						   _LU.getPointer(),_LU.getStride(), 
						   &_P[0], FFLAPACK::FflapackLQUP, &_Q[0] );
			
		}

		// Contruction of LQUP factorization of A (in-place in A)
		LQUPMatrix (const Field& F, BlasMatrix<Element>& A)
			: _F(F), _LU(A) , P(A.coldim()), Q(A.rowdim()), 
			  _m(A.rowdim()), _n(A.coldim())  {

			_rank= FFLAPACK::LUdivine( _F,FFLAS::FflasNonUnit, m, n, 
						   _LU.getPointer(),_LU.getStride(), 
						   &_P[0], FFLAPACK::FflapackLQUP, &_Q[0] );
			
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
		const TriangularBlasMatrix<Element>& getL() const;

		// get the matrix U
		const TriangularBlasMatrix<Element>& getU() const;

		// get the matrix S (from the LSP factorization of A deduced from LQUP)
		const BlasMatrix<Element>& getS() const;

		/*
		 * Solvers with Matrices
		 * Operand can be a BlasMatrix<Element> or a std::vector<Element>
		 */
		// solve AX=B
		template <class Operand>
		bool left_solve(Operand& X, const Operand& B) const {
			
			if ( ( _m != _n ) || ( _rank != _n ) || ( B.rowdim() != _n) )
				return false;
			X = *(new Operand(B) );
			return FactorizedMatrixLeftSolve<Field,Operand>()( _F, _LU, X );
		}

		// solve AX=B (X is stored in B)
		template <class Operand>
		bool left_solve(Operand& B) const {
			
			if ( ( _m != _n ) || ( _rank != _n ) || ( B.rowdim() != _n) )
				return false;
			return FactorizedMatrixLeftSolve<Field,Operand>()( _F, _LU, B );
		
		}

		// solve XA=B
		template <class Operand>
		bool right_solve(Operand& X, const Operand& B) const {
		
			if ( ( _m != _n ) || ( _rank != _n ) || ( B.coldim() != _m) )
				return false;
			X = *(new Operand(B) );
			return FactorizedMatrixRightSolve<Field,Operand>()( _F, _LU, X );
		}
		
		// solve XA=B (X is stored in B)
		template <class Operand>
		bool right_solve(Operand& B) const{

			if ( ( _m != _n ) || ( _rank != _n ) || ( B.coldim() != _m) )
				return false;
			return FactorizedMatrixRightSolve<Field,Operand>()( _F, _LU, B );
		}

		
		// solve LX=B (L from LQUP)
		template <class Operand>
		bool left_Lsolve(Operand& X, const Operand& B) const{
			if ( B.rowdim() != _m )
				return false;
			X = *(new Operand(B) );
			return FactorizedMatrixLeftLSolve<Field,Operand>()( _F, _LU, X, _rank );
		}
		
		// solve LX=B (L from LQUP) (X is stored in B)
		template <class Operand>
		bool left_Lsolve(Operand& B) const{
			if ( B.rowdim() != _m )
				return false;
			return FactorizedMatrixLeftLSolve<Field,Operand>()( _F, _LU, B, _rank );
		}

		// solve XL=B (L from LQUP)
		template <class Operand>
		bool right_Lsolve(Operand& X, const Operand& B) const{
			if ( B.coldim() != _m ) 
				return false;
			X = *(new Operand(B) );
			return FactorizedMatrixRightLSolve<Field,Operand>()( _F, _LU, X, _rank );
		}
		
		// solve XL=B (L from LQUP) (X is stored in B)
		template <class Operand>
		bool right_Lsolve(Operand& B) const{
			if ( B.coldim() != _m )
				return false;
			return FactorizedMatrixRightLSolve<Field,Operand>()( _F, _LU, B, _rank );
		}
		
		// solve UX=B (U from LQUP is r by r)
		template <class Operand>
		bool left_Usolve(Operand& X, const Operand& B) const{
			if ( B.rowdim() != _rank )
				return false;
			X = *(new Operand(B) );
			return FactorizedMatrixLeftUSolve<Field,Operand>()( _F, _LU, X );
		}
		
		// solve UX=B (U from LQUP) (X is stored in B)
		template <class Operand>
		bool rleft_Usolve(Operand& B) const{
			if ( B.rowdim() != _rank )
				return false;
			return FactorizedMatrixLeftUSolve<Field,Operand>()( _F, _LU, B );
		}

		// solve XU=B (U from LQUP)
		template <class Operand>
		bool right_Usolve(Operand& X, const Operand& B) const{
			if ( B.coldim() != _rank )
				return false;
			X = *(new Operand(B) );
			return FactorizedMatrixRightUSolve<Field,Operand>()( _F, _LU, X );
		}
		
		// solve XU=B (U from LQUP) (X is stored in B)
		template <class Operand>
		bool right_Usolve(Operand& B) const{
			if ( B.coldim() != _rank )
				return false;
			return FactorizedMatrixRightUSolve<Field,Operand>()( _F, _LU, B );
		}


	}; // end of class LQUPMatrix

} //end of namespace LinBox

#include <linbox/matrix/matrix-factorized.inl>

#endif
