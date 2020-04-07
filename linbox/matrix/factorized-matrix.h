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

#include <fflas-ffpack/fflas-ffpack.h>

#include "linbox/matrix/densematrix/blas-matrix.h"

#include "linbox/matrix/permutation-matrix.h"

namespace LinBox
{

	/*- @name Factorized Matrix
	 * @brief Solving using blas and LU style factored matrix.
	 */
	//-{

	// forward definition
	template <class Field>
	class PLUQMatrix;


	/*! PLUQ factorisation.
	 * This is a class to ease the use LU factorisation (see FFPACK::PLUQ
	 *
	 * The factorisation is \f$ A = P L U P \f$ with \c L lower unit
	 * triangular, \c U upper non-unit triangular, \c P and \c Q
	 * permutations.
	 *
	 * There are two kind of contructors (with and without permutations)
	 * and they build a \c PLUQ factorisation of a \c BlasMatrix/\c BlasBlackbox on
	 * a finite field.  There are methods for retrieving \p P \p L,\p U and \p Q
	 * matrices and methods for solving systems.
	 */
	//! @bug Should really be tempalted by Matrix and be a (sub)domain
	template <class Field>
	class PLUQMatrix {

	public:
		typedef typename Field::Element Element;
		//typedef std::vector<size_t> BlasPermutation;
        typedef  BlasMatrix<Field,typename Vector<Field>::Dense > matrixType;
        typedef  TriangularBlasMatrix<matrixType>  TriangularMatrix;
	protected:

		Field                     _field;
		BlasMatrix<Field,typename Vector<Field>::Dense >       &_factLU;
		BlasPermutation<size_t> &_permP;   //note: this is actually P^T!
		BlasPermutation<size_t> &_permQ;
		size_t                    _m;
		size_t                    _n;
		size_t                 _rank;
		bool                  _alloc;
		bool                  _plloc;

	public:

		//! Contruction of PLUQ factorization of A (making a copy of A)
		template<class _Rep>
		PLUQMatrix (const BlasMatrix<Field,_Rep>& A) ;

		//! Contruction of PLUQ factorization of A (in-place in A) 
		template<class _Rep>
		PLUQMatrix (BlasMatrix<Field,_Rep>& A) ;


		/*! Contruction of PLUQ factorization of A (making a copy of A).
		 * P and Q are arguments !
		 */
		template<class _Rep>
		PLUQMatrix (const BlasMatrix<Field,_Rep>& A,
			    BlasPermutation<size_t> & P, BlasPermutation<size_t> & Q) ;

		/*! Contruction of PLUQ factorization of A (in-place in A).
		 * P and Q are arguments !
		 * @bug in place ?
		 */
		template<class _Rep>
		PLUQMatrix (BlasMatrix<Field,_Rep>& A,
			    BlasPermutation<size_t> & P, BlasPermutation<size_t> & Q) ;

		//! destructor.
		~PLUQMatrix () ;

		//! get the field on which the factorization is done
		Field& field() ;

		//! get row dimension
		size_t rowdim() const ;

		//! get column dimension
		size_t coldim() const ;

		//! get the rank of matrix
		size_t getRank() const ;

		/*! get the permutation Q.
		 * (no copy)
		 */
		const BlasPermutation<size_t>& getQ() const ;

		/*! get the permutation Q.
		 * (copy)
		 */
		BlasPermutation<size_t> & getQ( BlasPermutation<size_t> & Q ) const ;

		/** Get the <i>transpose</i> of the permutation \p P.
		 * @warning This does not return \p P itself! (because it is
		 * more difficult to compute) If needed, \p P can be obtained
		 * as a \p TransposedBlasMatrix from the return value. One
		 * reason this confusion exists is that left-multiplying by
		 * a permuation matrix corresponds to a row permuation \f$\pi \in S_n\f$,
		 * while right-multiplying by the same matrix corresponds to
		 * the inverse column permutation \f$\pi^{-1}\f$!  Usually this
		 * is handled intelligently (eg by \c applyP) but you must be
		 * careful with \c getP().
		 */
		const BlasPermutation<size_t>& getP() const ;

		/*! get the permutation P^T.
		 * (copy)
		 * @warning see <code>PLUQMatrix::getP()</code>
		 */
		BlasPermutation<size_t> & getP( BlasPermutation<size_t> & PT ) const ;

		/*! get the Matrix \c  L.
		 * @param[out] L
		 * @param _QLUP if true then \c L form \c QLUP decomposition,
		 * else \c L is form \c PLUQ decomposition.
		 * @pre \c L has unit diagonal
		 */
        
        TriangularMatrix& getL(TriangularMatrix& L, bool _QLUP = false) const;

		/*! get the matrix  \c  U.
		 * @pre \c   U has non-unit diagonal
		 */
		TriangularMatrix& getU(TriangularMatrix& U) const;

		// /*! get the matrix S.
		//  *  from the LSP factorization of A deduced from PLUQ)
		//  */
		// template<class _Rep>
		// BlasMatrix<Field,_Rep>& getS( BlasMatrix<Field,_Rep>& S) const;

		/*! @internal get a pointer to the begin of storage.
		*/
		Element* getPointer() const ;

		/*! @internal get  the stride in \c _factLU
		*/
		size_t getStride() const ;

		/*!
		 * Solvers with matrices or vectors
		 * Operand can be a BlasMatrix<Field,_Rep> or a std::vector<Element>
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

		// solve LX=B (L from PLUQ)
		template <class Operand>
		Operand& left_Lsolve(Operand& X, const Operand& B) const;

		// solve LX=B (L from PLUQ) (X is stored in B)
		template <class Operand>
		Operand& left_Lsolve(Operand& B) const;

		// solve XL=B (L from PLUQ)
		template <class Operand>
		Operand& right_Lsolve(Operand& X, const Operand& B) const;

		// solve XL=B (L from PLUQ) (X is stored in B)
		template <class Operand>
		Operand& right_Lsolve(Operand& B) const;

		// solve UX=B (U from PLUQ is r by r)
		template <class Operand>
		Operand& left_Usolve(Operand& X, const Operand& B) const;

		// solve UX=B (U from PLUQ) (X is stored in B)
		template <class Operand>
		Operand& rleft_Usolve(Operand& B) const;

		// solve XU=B (U from PLUQ)
		template <class Operand>
		Operand& right_Usolve(Operand& X, const Operand& B) const;

		// solve XU=B (U from PLUQ) (X is stored in B)
		template <class Operand>
		Operand& right_Usolve(Operand& B) const;
		//@}


	}; // end of class PLUQMatrix

	//-}
} // end of namespace LinBox

#include "linbox/matrix/factorized-matrix.inl"

#endif //__LINBOX_factorized_matrix_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
