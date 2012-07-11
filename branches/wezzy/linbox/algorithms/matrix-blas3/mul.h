/*  Copyright (C) 2012 the members of the LinBox group
 * Written by B. Boyer < bboyer@imag.fr >
 *
 * This file is part of the LinBox library.
 *
 * ========LICENCE========
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * LinBox is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 */

/*! @file matrix-blas3/mul.h
 * @ingroup algorithm
 * @ingroup blas
 * @brief BLAS3 multiplication algorithms.
 */

#ifndef __LINBOX_matrix_blas3_mul_H
#define __LINBOX_matrix_blas3_mul_H

#include <vector>
#include "linbox/integer.h"
#include <linbox/field/givaro.h>
#include <linbox/matrix/blas-matrix.h>



// Methods
namespace LinBox {
	namespace BLAS3 {

		//! BLAS3 Multiplication methods
		namespace mulMethod {
			template<class Field>
			struct ToomCook {
				const Field & _myF;
				bool memory_unlimited ;
				ToomCook(Field & myF, bool mem = false) :
					_myF(myF),
					memory_unlimited(mem)
				{};
			// this is the place for the helper/caching (TC,iTC).
			} ; //! Toom-Cook method.

			struct naive {};

			struct FLINT {};
		}
	}
}

// ToomCook
namespace LinBox {
	namespace BLAS3 {

		/** @brief Build the Toom-Cook matrix helper.
		 * We only provide this function for dense matrices (Blas for
		 * the moment) over Z/pZ.
		 * @param [out] TC matrix of Fp(i^j)
		 * @param [out] iTC inverse of TC
		 * @return \p TC
		 */
		template<class Zpz>
		BlasMatrix<Zpz>& ToomCook(BlasMatrix<Zpz>& TC, BlasMatrix<Zpz>& iTC);

		/** @brief Toom-Cook multiplication for GF(p^e)
		 * A matrix over GF(p^e) is represented by a polynomial of Z/pZ matrices.
		 * @tparam Zpz is some modular field, GFq is GivaroExtension (or the like)
		 * @param [out] C result
		 * @param A matrix
		 * @param B matrix
		 * @return C=AB
                 * @warning p should not be too small, and e>1 (you've been warned...)
		 */
		template<class Zpz, class GFq>
		std::vector<BlasMatrix<Zpz> >& mul (std::vector<BlasMatrix<Zpz> >& C,
						    const std::vector<BlasMatrix<Zpz> >&A,
						    const std::vector<BlasMatrix<Zpz> >&B,
						    const mulMethod::ToomCook<GFq>& T);

		/** @brief Toom-Cook multiplication for GF(p^e)
                 *
		 * @tparam Zpz is some modular field, GFq is GivaroExtension (or the like)
		 * @param [out] C result
		 * @param A matrix
		 * @param B matrix
		 * @return C=AB
                 * @warning p should not be too small, and e>1 (you've been warned...)
		 */
		template<class Zpz>
		BlasMatrix<GivaroExtension<Zpz> >&
		mul (BlasMatrix<GivaroExtension<Zpz> >& C,
			 const BlasMatrix<GivaroExtension<Zpz> >& A,
			 const BlasMatrix<GivaroExtension<Zpz> >& B,
			 const mulMethod::ToomCook<GivaroExtension<Zpz> > & T);

#if 0 /* Generic method */
		template<class ZpzMatrix>
		std::vector<ZpzMatrix >& mul (std::vector<ZpzMatrix >& C,
					      std::vector<ZpzMatrix >&A,
					      std::vector<ZpzMatrix >&B,
					      const mulMethod::ToomCook &);
#endif



	} // BLAS3
}


#include "linbox/algorithms/matrix-blas3/mul-toomcook.inl"

// naive
namespace LinBox {
	namespace BLAS3 {

		/** @brief Triple loop !
		 * just a simple triple loop
		 * @pre works if _anyMatrix has setEntry and _otherMatrix has getEntry
		 * @param [out] C result
		 * @param [in] A matrix
		 * @param [in] B matrix
		 * @return C=AB
		 */
		template<class _anyMatrix, class _otherMatrix1, class _otherMatrix2>
		_anyMatrix & mul (_anyMatrix& C,
				  const _otherMatrix1& A,
				  const _otherMatrix2& B,
				  const mulMethod::naive &);


	}
}
#include "linbox/algorithms/matrix-blas3/mul-naive.inl"

// flint

#ifdef __LINBOX_HAVE_FLINT
namespace FLINT {
	extern "C" {
		// #include "longlong.h"
#define __GMP_BITS_PER_MP_LIMB GMP_LIMB_BITS
#include "flint.h"
#include "fmpz_mat.h"
	}
}
namespace LinBox {
	namespace BLAS3 {

		//! Wrapper to FLINT mat-mul
		template<class DenseIntMat>
		DenseIntMat &
		mul (DenseIntMat& C,
			 const DenseIntMat& A,
			 const DenseIntMat& B,
			 const mulMethod::FLINT & );

	}
}
#endif // __LINBOX_HAVE_FLINT
//#include "linbox/algorithms/matrix-blas3/mul-flint.inl"



// <+other algo+>
namespace LinBox {
	namespace BLAS3 {

	}
}
//#include "linbox/algorithms/matrix-blas3/mul-<+other algo+>.inl"

#endif // __LINBOX_matrix_blas3_mul_H

//Local Variables:
//mode: C++
//tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

