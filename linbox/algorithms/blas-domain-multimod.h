/* linbox/algorithms/blas-domain.inl
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



#ifndef __LINBOX_blas_matrix_domain_multimod_H
#define __LINBOX_blas_matrix_domain_multimod_H

#include "linbox/matrix/blas-matrix-multimod.h"
#include "linbox/algorithms/blas-domain.h"

THIS_CODE_MAY_NOT_COMPILE_AND_IS_NOT_TESTED

namespace LinBox
{


	/*! specialisation for MultiModDouble.
	*/
// #ifndef __INTEL_COMPILER
	// template <>
// #endif
	class BlasMatrixDomainInv<MultiModDouble,BlasMatrix<MultiModDouble> > {
	public:
		int operator() (const MultiModDouble                   &F,
				BlasMatrix<MultiModDouble>        &Ainv,
				const BlasMatrix<MultiModDouble>     &A) const
		{

			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.rowdim() == Ainv.rowdim());
			linbox_check( A.coldim() == Ainv.coldim());
			BlasMatrix<MultiModDouble> tmp(A);
			return (*this)(F,Ainv,tmp);
		}

		int operator() (const MultiModDouble                &F,
				BlasMatrix<MultiModDouble>     &Ainv,
				BlasMatrix<MultiModDouble>        &A) const
		{

			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.rowdim() == Ainv.rowdim());
			linbox_check( A.coldim() == Ainv.coldim());
			int nullity, defrank=0;

			for (size_t i=0;i<F.size();++i){
				FFPACK::Invert(F.getBase(i),A.rowdim(), A.getMatrix(i)->getPointer(),A.getMatrix(i)->getStride(),
					       Ainv.getMatrix(i)->getPointer(),Ainv.getMatrix(i)->getStride(),nullity);
				defrank+=nullity;
			}
			return defrank;
		}

	};

}


#endif // __LINBOX_blas_matrix_domain_multimod_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

