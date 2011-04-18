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


#ifndef __LINBOX_factorized_matrix_H
#define __LINBOX_factorized_matrix_H

#if 0
#undef _P
#undef _Q
#endif

#include <vector>

#include <linbox/matrix/blas-matrix.h>
#include <linbox/blackbox/blas-blackbox.h>
#include <linbox/algorithms/blas-domain.h>
#include <linbox/ffpack/ffpack.h>

namespace LinBox
{

/** @name Factorized Matrix
 * @brief Solving using blas and LU style factored matrix.
 */
//@{

	template <class Field>
	class LQUPMatrix;

	/// 
	template <class Field, class Operand> 
	class FactorizedMatrixLeftSolve {
	public:
		Operand& operator() ( const Field& F, 
				      const LQUPMatrix<Field>& A, 
				      Operand& X, const Operand& B ) const;
		Operand& operator() ( const Field& F,
				      const LQUPMatrix<Field>& A, 
				      Operand& B ) const;
	}; // end of class FactorizedMatrixLeftSolve

	/// 
	template <class Field, class Operand> 
	class FactorizedMatrixRightSolve {
	public:
		Operand& operator() ( const Field& F, 
				      const LQUPMatrix<Field>& A, 
				      Operand& X, const Operand& B ) const;
		Operand& operator() ( const Field& F, 
				      const LQUPMatrix<Field>& A, 
				      Operand& B ) const;
	}; // end of class FactorizedMatrixRightSolve

	///
	template <class Field, class Operand> 
	class FactorizedMatrixLeftLSolve {
	public:
		Operand& operator() ( const Field& F, 
				      const LQUPMatrix<Field>& A,
				      Operand& X, const Operand& B ) const;
		Operand& operator() ( const Field& F, 
				      const LQUPMatrix<Field>& A,
				      Operand& B ) const;
	}; // end of class FactorizedMatrixLeftLSolve

	///
	template <class Field, class Operand> 
	class FactorizedMatrixRightLSolve {
	public:
		Operand& operator() ( const Field& F, 
				      const LQUPMatrix<Field>& A,
				      Operand& X, const Operand& B ) const;
		Operand& operator() ( const Field& F, 
				      const LQUPMatrix<Field>& A,
				      Operand& B ) const;
	}; // end of class FactorizedMatrixRightLsolve

	///
	template <class Field, class Operand> 
	class FactorizedMatrixLeftUSolve {
	public:
		Operand& operator() ( const Field& F, 
				      const LQUPMatrix<Field>& A,
				      Operand& X, const Operand& B ) const;
		Operand& operator() ( const Field& F, 
				      const LQUPMatrix<Field>& A,
				      Operand& B ) const;
	}; // end of class FactorizedMatrixLeftUSolve

	///
	template <class Field, class Operand> 
	class FactorizedMatrixRightUSolve {
	public:
		Operand& operator() ( const Field& F, 
				      const LQUPMatrix<Field>& A,
				      Operand& X, const Operand& B ) const;
		Operand& operator() ( const Field& F, 
				      const LQUPMatrix<Field>& A,
				      Operand& B ) const;
	}; // end of class FactorizedMatrixRightUSolve
	
	///
	template <class Field>
	class LQUPMatrix {

	public:
		typedef typename Field::Element Element;
		//typedef std::vector<size_t> BlasPermutation;

	protected:

		Field                   _F;
		BlasMatrix<Element>    &_LU;
		BlasPermutation         _PP;
		BlasPermutation         _QQ;  //note: this is actually Qt! 
		size_t                  _m;
		size_t                  _n;
		size_t               _rank;
		bool                 _alloc;
    
	public:


		// Contruction of LQUP factorization of A (making a copy of A)
		LQUPMatrix (const Field& F, const BlasMatrix<Element>& A)
			: _F(F), _LU(*(new BlasMatrix<Element> (A))) ,
			  _PP(A.coldim()), _QQ(A.rowdim()), _m(A.rowdim()),
			  _n(A.coldim()), _alloc(true)  {
			//std::cerr<<"Je passe par le constructeur const"<<std::endl;

			_rank= FFPACK::LUdivine( _F,FFLAS::FflasNonUnit,  FFLAS::FflasNoTrans, _m, _n, 
						 _LU.getPointer(),_LU.getStride(), 
						 _PP.getWritePointer(), _QQ.getWritePointer(), FFPACK::FfpackLQUP );
			
		}

		// Contruction of LQUP factorization of A (in-place in A)
		LQUPMatrix (const Field& F, BlasMatrix<Element>& A)
			: _F(F), _LU(A) , _PP(A.coldim()), _QQ(A.rowdim()), 
			  _m(A.rowdim()), _n(A.coldim()), _alloc(false) {
			//std::cerr<<"Je passe par le constructeur non const"<<std::endl;
			_rank= FFPACK::LUdivine( _F,FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, _m, _n, 
						 _LU.getPointer(),_LU.getStride(), 
						 _PP.getWritePointer(), _QQ.getWritePointer(), FFPACK::FfpackLQUP );
			
		}

		// Contruction of LQUP factorization of A (making a copy of A)
		LQUPMatrix (const BlasBlackbox<Field>& A)
			: _F(A.field()), _LU(*(new BlasMatrix<Element> (A))) ,
			  _PP(A.coldim()), _QQ(A.rowdim()), _m(A.rowdim()),
			  _n(A.coldim()), _alloc(true)  {
			
			_rank= FFPACK::LUdivine( _F,FFLAS::FflasNonUnit,  FFLAS::FflasNoTrans,  _m, _n, 
						 _LU.getPointer(),_LU.getStride(), 
						 _PP.getWritePointer(), _QQ.getWritePointer(), FFPACK::FfpackLQUP );
		}

		// Contruction of LQUP factorization of A (in-place in A)
		LQUPMatrix (BlasBlackbox<Field>& A)
			: _F(A.field()), _LU(static_cast<BlasMatrix<Element>&> (A)) , _PP(A.coldim()), _QQ(A.rowdim()), 
			  _m(A.rowdim()), _n(A.coldim()), _alloc(false) {
			
			_rank= FFPACK::LUdivine( _F,FFLAS::FflasNonUnit,  FFLAS::FflasNoTrans, _m, _n, 
						 _LU.getPointer(),_LU.getStride(), 
						 _PP.getWritePointer(), _QQ.getWritePointer(), FFPACK::FfpackLQUP );
		}
		

		~LQUPMatrix () {
			if (_alloc)
				delete &_LU;
		}

		// get the field on which the factorization is done
		Field& field() {return _F;}

		// get row dimension
		size_t rowdim() const {return _m;}

		// get column dimension
		size_t coldim() const {return _n;}
    
		// get the rank of matrix
		size_t getrank() const {return _rank;}
    
		// get the permutation P
		const BlasPermutation& getP() const {return _PP;}
       
		/** get the _transpose_ of the permutation Q
		 * NOTE: this does not return Q itself! (because it is more difficult to compute)
		 * If needed, Q can be obtained as a TransposedBlasMatrix from the return value
		 *
		 * One reason this confusion exists is that left-multiplying by a permuation matrix 
		 * corresponds to a row permuation \pi \in S_n, while right-multiplying by the same matrix
		 * corresponds to the inverse column permutation \pi^(-1)!
		 * Usually this is handled intelligently (eg by applyP) but you must be careful with getQ().
		 */
		const BlasPermutation& getQ() const  {return _QQ;}

		// get the Matrix L
		TriangularBlasMatrix<Element>& getL(TriangularBlasMatrix<Element>& L) const;

		// get the matrix U
		TriangularBlasMatrix<Element>& getU(TriangularBlasMatrix<Element>& U) const;

		// get the matrix S (from the LSP factorization of A deduced from LQUP)
		BlasMatrix<Element>& getS( BlasMatrix<Element>& S) const;

		// get a pointer to the begin of storage
		Element* getPointer() const { return _LU.getPointer(); }

		// get a pointer to the begin of storage
		//const size_t getStride() const { return _LU.getStride(); }
		size_t getStride() const { return _LU.getStride(); } //BB: pas besoin

		/*
		 * Solvers with matrices or vectors
		 * Operand can be a BlasMatrix<Element> or a std::vector<Element>
		 */
		// solve AX=B
		template <class Operand>
		Operand& left_solve(Operand& X, const Operand& B) const {
			return FactorizedMatrixLeftSolve<Field,Operand>()( _F, *this, X, B );
		}

		// solve AX=B (X is stored in B)
		template <class Operand>
		Operand& left_solve(Operand& B) const {
			return FactorizedMatrixLeftSolve<Field,Operand>()( _F, *this, B );
		}

		// solve XA=B
		template <class Operand>
		Operand& right_solve(Operand& X, const Operand& B) const {
			return FactorizedMatrixRightSolve<Field,Operand>()( _F, *this, X, B );
		}
		
		// solve XA=B (X is stored in B)
		template <class Operand>
		Operand& right_solve(Operand& B) const{
			return FactorizedMatrixRightSolve<Field,Operand>()( _F, *this, B );
		}
		
		// solve LX=B (L from LQUP)
		template <class Operand>
		Operand& left_Lsolve(Operand& X, const Operand& B) const{
			return FactorizedMatrixLeftLSolve<Field,Operand>()( _F, *this, X, B );
		}
		
		// solve LX=B (L from LQUP) (X is stored in B)
		template <class Operand>
		Operand& left_Lsolve(Operand& B) const{
			return FactorizedMatrixLeftLSolve<Field,Operand>()( _F, *this, B );
		}

		// solve XL=B (L from LQUP)
		template <class Operand>
		Operand& right_Lsolve(Operand& X, const Operand& B) const{
			return FactorizedMatrixRightLSolve<Field,Operand>()( _F, *this, X, B );
		}
		
		// solve XL=B (L from LQUP) (X is stored in B)
		template <class Operand>
		Operand& right_Lsolve(Operand& B) const{
			return FactorizedMatrixRightLSolve<Field,Operand>()( _F, *this, B );
		}
		
		// solve UX=B (U from LQUP is r by r)
		template <class Operand>
		Operand& left_Usolve(Operand& X, const Operand& B) const{
			return FactorizedMatrixLeftUSolve<Field,Operand>()( _F, *this, X, B );
		}
		
		// solve UX=B (U from LQUP) (X is stored in B)
		template <class Operand>
		Operand& rleft_Usolve(Operand& B) const{
			return FactorizedMatrixLeftUSolve<Field,Operand>()( _F, *this, B );
		}

		// solve XU=B (U from LQUP)
		template <class Operand>
		Operand& right_Usolve(Operand& X, const Operand& B) const{
			return FactorizedMatrixRightUSolve<Field,Operand>()( _F, *this, X, B );
		}
		
		// solve XU=B (U from LQUP) (X is stored in B)
		template <class Operand>
		Operand& right_Usolve(Operand& B) const{
			return FactorizedMatrixRightUSolve<Field,Operand>()( _F, *this, B );
		}


	}; // end of class LQUPMatrix

//@}
} // end of namespace LinBox

#include <linbox/matrix/factorized-matrix.inl>

#endif //__LINBOX_factorized_matrix_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
