/* Copyright (C) 2010 LinBox
 * Written by Brice Boyer <brice.boyer@imag.fr>
 *
 *
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

#ifndef __LINBOX_matrix_random_matrix_H
#define __LINBOX_matrix_random_matrix_H

/** @file matrix/random-matrix.h
 * @ingroup matrix
 * @brief Implementation of random matrices.
 *
 * We provide function to create random matrices (dense, sparse, structured)
 * on several rings. This header was first introduced to avoid code redundancy in tests/
 * and make it easier to write tests/ examples/.
 *
 * @todo à la vector/stream.h
 * @bug this belongs to algorithms...
 */

#include "linbox/matrix/dense-matrix.h"
#include "linbox/randiter/random-integer.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/matrix/permutation-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/cra-full-multip-fixed.h"

namespace LinBox
{
	template<class Field> class BlasMatrixDomain ;

	// non dependant struct -> out of class.
	struct RankBuildMethod {
			RankBuildMethod() {}
		} ;
		struct	_LU_ : public RankBuildMethod {
			_LU_(){}
		} ;
		// LU_sparse_,
		// LU_cra_,
		struct _Rank_update_ : public RankBuildMethod {
			_Rank_update_() {}
		} ;

		/// random method for constructing rank
		struct RankBuilder {
			// private :
			// balancedLC_
			// public:
			typedef _LU_                   LU_ ;
			typedef _Rank_update_ Rank_update_ ;
			RankBuilder(){}
		};

	/// Random Dense Matrix builder.
	template<class Randiter, class Field>
	class RandomDenseMatrix {


	protected:

	private :
		Field      F_ ; //!< The field containing the random entries. @todo is there a copy made ?
		/*! How are entries generated ?
		 * @pre need only provide <code>elmt& random(elmt&);</code>
		 * @see \ref LinBox::RandIterArchetype
		 */
		Randiter   R_ ;
		// Matrix    & A_ ; //!< The resulting random matrix
	public :
		/// constructor
		RandomDenseMatrix(Field & F, Randiter & R) :
			F_(F), R_(R)
		{  }

		RandomDenseMatrix(const Field & F, Randiter & R) :
			F_(F), R_(R)
		{  }

		/// destructor
		~RandomDenseMatrix() {}

		/*! creates a randomly filled matrix.
		 * @param A matrix to be randomized.
		 */
		template<class Matrix>
		Matrix & random(Matrix & A) ;

		/*! provide a matrix with prescribed rank.
		 * Default method.
		 * @param A
		 * @param rank expected rank
		 * @warning No certificate yet.
		 */
		template<class Matrix>
		Matrix & randomRank(Matrix & A, int rank);


		/*! provide a matrix with prescribed rank.
		 * @param A
		 * @param rank expected rank
		 * @param meth how is the matrix generated ? see \ref RankBuilder.
		 * @warning No certificate yet.
		 */

		template<class Matrix>
		Matrix & randomRank(Matrix & A, int rank
				    , const RankBuilder::LU_ & meth );

		template<class Matrix>
		Matrix & randomRank(Matrix & A, int rank
				    , const RankBuilder::Rank_update_ & meth );



		// template<class Matrix>
		// void randomInvertible();

		// void randomNilpotent(int nil k); // P N_k P^(-1)
		//
		// Matrix& getMatrix(Matrix & A) {
		// return A = A_ ;
		// }



	};

	/// @todo To be factorized.
	void RandomBlasPermutation(BlasPermutation<size_t> & P) ;
} // LinBox

#include "linbox/matrix/random-matrix.inl"

#endif // __LINBOX_matrix_random_matrix_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

