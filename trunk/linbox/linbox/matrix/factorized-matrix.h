/* linbox/matrix/factorized-matrix.h
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@ens-lyon.fr
 *               Clément Pernet clement.pernet@imag.fr
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


#ifndef __LINBOX_factorized_matrix_H
#define __LINBOX_factorized_matrix_H


#include <vector>

#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blas-domain.h"
#include <fflas-ffpack/ffpack/ffpack.h>

#include "linbox/matrix/permutation-matrix.h"

namespace LinBox
{

	/*- @name Factorized Matrix
	 * @brief Solving using blas and LU style factored matrix.
	 */
	//-{

	// forward definition
	template <class Field>
	class LQUPMatrix;


	/*! LQUP factorisation.
	 * This is a class to ease the use LU factorisation (see \ref LUdivine
	 * in \ref FFPACK)
	 *
	 * The factorisation is \f$ A = L Q U P \f$ with \c L lower unit
	 * triangular, \c U upper non-unit triangular, \c P and \c Q
	 * permutations.
	 *
	 * There are two kind of contructors (with and without permutations)
	 * and they build a \c LQUP factorisation of a \c BlasMatrix/\c BlasBlackbox on
	 * a finite field.  There are methods for retrieving \p L,\p Q,\p U and \p P
	 * matrices and methods for solving systems.
	 */
	template <class Field>
	class LQUPMatrix {

	public:
		typedef typename Field::Element Element;
		//typedef std::vector<size_t> BlasPermutation;

	protected:

		Field                     _field;
		BlasMatrix<Field>       &_factLU;
		BlasPermutation<size_t> &_permP;
		BlasPermutation<size_t> &_permQ;  //note: this is actually Qt!
		size_t                    _m;
		size_t                    _n;
		size_t                 _rank;
		bool                  _alloc;
		bool                  _plloc;

	public:

#if 0
		//! Contruction of LQUP factorization of A (making a copy of A)
		LQUPMatrix (const Field& F, const BlasMatrix<Field>& A) :
			_field(F), _factLU(*(new BlasMatrix<Field> (A))) ,
			_permP(*(new BlasPermutation<size_t>(A.coldim()))),
			_permQ(*(new BlasPermutation<size_t>(A.rowdim()))),
			_m(A.rowdim()), _n(A.coldim()),
			_alloc(true),_plloc(true)
		{
			//std::cerr<<"Je passe par le constructeur const"<<std::endl;

			_rank= FFPACK::LUdivine( _field,FFLAS::FflasNonUnit,  FFLAS::FflasNoTrans, _m, _n,
						 _factLU.getPointer(),_factLU.getStride(),
						 _permP.getWritePointer(), _permQ.getWritePointer(), FFPACK::FfpackLQUP );
			_permP.setOrder(_rank);
			_permQ.setOrder(_rank);

		}

		//! Contruction of LQUP factorization of A (in-place in A)
		LQUPMatrix (const Field& F, BlasMatrix<Field>& A) :
			_field(F), _factLU(A) ,
			_permP(*(new BlasPermutation<size_t>(A.coldim()))),
			_permQ(*(new BlasPermutation<size_t>(A.rowdim()))),
			_m(A.rowdim()), _n(A.coldim()),
			_alloc(false),_plloc(true)
		{
			if (!A.coldim() || !A.rowdim()) {
				// throw LinBoxError("LQUP does not accept empty matrices");
				_rank = 0 ;
			}
			else {
				//std::cerr<<"Je passe par le constructeur non const"<<std::endl;
				_rank= FFPACK::LUdivine( _field,FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, _m, _n,
							 _factLU.getPointer(),_factLU.getStride(),
							 _permP.getWritePointer(), _permQ.getWritePointer(), FFPACK::FfpackLQUP );
			}
			_permP.setOrder(_rank);
			_permQ.setOrder(_rank);

		}

		/*! Contruction of LQUP factorization of A (making a copy of A).
		 * P and Q are arguments !
		 */
		LQUPMatrix (const BlasMatrix<Field>& A,
			    BlasPermutation<size_t> & P, BlasPermutation<size_t> & Q) :
			_field(F), _factLU(*(new BlasMatrix<Field> (A))) ,
			_permP(P), _permQ(Q),
			_m(A.rowdim()), _n(A.coldim()),
			_alloc(true),_plloc(false)
		{
			//std::cerr<<"Je passe par le constructeur const"<<std::endl;

			linbox_check(_permQ.getOrder()==A.rowdim());
			linbox_check(_permP.getOrder()==A.coldim());
			_rank= FFPACK::LUdivine( _field,FFLAS::FflasNonUnit,  FFLAS::FflasNoTrans, _m, _n,
						 _factLU.getPointer(),_factLU.getStride(),
						 _permP.getWritePointer(), _permQ.getWritePointer(), FFPACK::FfpackLQUP );
			_permP.setOrder(_rank);
			_permQ.setOrder(_rank);


		}

		/*! Contruction of LQUP factorization of A (in-place in A).
		 * P and Q are arguments !
		 */
		LQUPMatrix ( BlasMatrix<Field>& A,
			    BlasPermutation<size_t> & P, BlasPermutation<size_t> & Q) :
			_field(F), _factLU(A) , _permP(P), _permQ(Q),
			_m(A.rowdim()), _n(A.coldim()),
			_alloc(false),_plloc(false)
		{
			//std::cerr<<"Je passe par le constructeur non const"<<std::endl;
			linbox_check(_permQ.getOrder()<=A.rowdim());
			linbox_check(_permP.getOrder()<=A.coldim());
			if (_permQ.getOrder() == 0)
				_permQ.resize(A.rowdim());
			if (_permP.getOrder() == 0)
				_permP.resize(A.coldim());

			_rank= FFPACK::LUdivine( _field,FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, _m, _n,
						 _factLU.getPointer(),_factLU.getStride(),
						 _permP.getWritePointer(), _permQ.getWritePointer(), FFPACK::FfpackLQUP );
			_permP.setOrder(_rank);
			_permQ.setOrder(_rank);

		}
#endif
		//! Contruction of LQUP factorization of A (making a copy of A)
		LQUPMatrix (const BlasMatrix<Field>& A) ;

		//! Contruction of LQUP factorization of A (in-place in A)
		LQUPMatrix (BlasMatrix<Field>& A) ;


		/*! Contruction of LQUP factorization of A (making a copy of A).
		 * P and Q are arguments !
		 */
		LQUPMatrix (const BlasMatrix<Field>& A,
			    BlasPermutation<size_t> & P, BlasPermutation<size_t> & Q) ;

		/*! Contruction of LQUP factorization of A (in-place in A).
		 * P and Q are arguments !
		 * @bug in place ?
		 */
		LQUPMatrix (BlasMatrix<Field>& A,
			    BlasPermutation<size_t> & P, BlasPermutation<size_t> & Q) ;

		//! destructor.
		~LQUPMatrix () ;

		//! get the field on which the factorization is done
		Field& field() ;

		//! get row dimension
		size_t rowdim() const ;

		//! get column dimension
		size_t coldim() const ;

		//! get the rank of matrix
		size_t getRank() const ;

		/*! get the permutation P.
		 * (no copy)
		 */
		const BlasPermutation<size_t>& getP() const ;

		/*! get the permutation P.
		 * (copy)
		 */
		BlasPermutation<size_t> & getP( BlasPermutation<size_t> & P ) const ;

		/** Get the <i>transpose</i> of the permutation \p Q.
		 * @warning This does not return \p Q itself! (because it is
		 * more difficult to compute) If needed, \p Q can be obtained
		 * as a \p TransposedBlasMatrix from the return value. One
		 * reason this confusion exists is that left-multiplying by
		 * a permuation matrix corresponds to a row permuation \f$\pi \in S_n\f$,
		 * while right-multiplying by the same matrix corresponds to
		 * the inverse column permutation \f$\pi^{-1}\f$!  Usually this
		 * is handled intelligently (eg by \c applyP) but you must be
		 * careful with \c getQ().
		 */
		const BlasPermutation<size_t>& getQ() const ;

		/*! get the permutation Qt.
		 * (copy)
		 * @warning see <code>LQUPMatrix::getQ()</code>
		 */
		BlasPermutation<size_t> & getQ( BlasPermutation<size_t> & Qt ) const ;

		/*! get the Matrix \c  L.
		 * @param[out] L
		 * @param _QLUP if true then \c L form \c QLUP decomposition,
		 * else \c L is form \c LQUP decomposition.
		 * @pre \c L has unit diagonal
		 */
		TriangularBlasMatrix<Field>& getL(TriangularBlasMatrix<Field>& L, bool _QLUP = false) const;

		/*! get the matrix  \c  U.
		 * @pre \c   U has non-unit diagonal
		 */
		TriangularBlasMatrix<Field>& getU(TriangularBlasMatrix<Field>& U) const;

		/*! get the matrix S.
		 *  from the LSP factorization of A deduced from LQUP)
		 */
		BlasMatrix<Field>& getS( BlasMatrix<Field>& S) const;

		/*! @internal get a pointer to the begin of storage.
		*/
		Element* getPointer() const ;

		/*! @internal get  the stride in \c _factLU
		*/
		size_t getStride() const ;

		/*!
		 * Solvers with matrices or vectors
		 * Operand can be a BlasMatrix<Field> or a std::vector<Element>
		 */
		//@{
		// solve AX=B
		template <class Operand>
		Operand& left_solve(Operand& X, const Operand& B) const;

		// solve AX=B (X is stored in B)
		template <class Operand>
		Operand& left_solve(Operand& B) const;

		// solve XA=B
		template <class Operand>
		Operand& right_solve(Operand& X, const Operand& B) const;

		// solve XA=B (X is stored in B)
		template <class Operand>
		Operand& right_solve(Operand& B) const;

		// solve LX=B (L from LQUP)
		template <class Operand>
		Operand& left_Lsolve(Operand& X, const Operand& B) const;

		// solve LX=B (L from LQUP) (X is stored in B)
		template <class Operand>
		Operand& left_Lsolve(Operand& B) const;

		// solve XL=B (L from LQUP)
		template <class Operand>
		Operand& right_Lsolve(Operand& X, const Operand& B) const;

		// solve XL=B (L from LQUP) (X is stored in B)
		template <class Operand>
		Operand& right_Lsolve(Operand& B) const;

		// solve UX=B (U from LQUP is r by r)
		template <class Operand>
		Operand& left_Usolve(Operand& X, const Operand& B) const;

		// solve UX=B (U from LQUP) (X is stored in B)
		template <class Operand>
		Operand& rleft_Usolve(Operand& B) const;

		// solve XU=B (U from LQUP)
		template <class Operand>
		Operand& right_Usolve(Operand& X, const Operand& B) const;

		// solve XU=B (U from LQUP) (X is stored in B)
		template <class Operand>
		Operand& right_Usolve(Operand& B) const;
		//@}


	}; // end of class LQUPMatrix

	//-}
} // end of namespace LinBox

#include "linbox/matrix/factorized-matrix.inl"

#endif //__LINBOX_factorized_matrix_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

