/*
 * Copyright (C) 2004,2014 Pascal Giorgi, Clément Pernet, the LinBox group
 *
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@ens-lyon.fr
 *               Clément Pernet clement.pernet@imag.fr
 *               Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

/*! @file matrix/densematrix/blas-tranposed-matrix.h
 * @ingroup densematrix
 * @brief Transposed Matrix
 * not much.
 *
 */

#ifndef __LINBOX_matrix_densematrix_blas_transposed_matrix_H
#define __LINBOX_matrix_densematrix_blas_transposed_matrix_H
#include "linbox/matrix/transpose-matrix.h"

namespace LinBox
{ /*  Transposed Matrix */    
	/*! TransposedBlasMatrix.
	 * NO DOC
	 */
	template< class Matrix >
	class TransposedBlasMatrix : public TransposeMatrix<Matrix> { // most methods are inherited 
	public:
        typedef TransposedBlasMatrix<Matrix>                Self_t;
        typedef typename Matrix::Field                       Field;
        typedef typename Field::Element                    Element;  
        typedef typename Field::Element_ptr            Element_ptr;      //!< Pointer to Element type
        typedef typename Field::ConstElement_ptr   ConstElement_ptr; //!< Pointer to const Element type
        
        using TransposeMatrix<Matrix>::TransposeMatrix; // use cstor from TransposeMatrix
        
        Self_t& copy(const Self_t& A) { return *this=A;}

        size_t getStride() const { return TransposeMatrix<Matrix>::_Mat.getStride();}
        
        ConstElement_ptr getPointer()      const { return TransposeMatrix<Matrix>::_Mat.getPointer();}
        Element_ptr      getPointer()            { return TransposeMatrix<Matrix>::_Mat.getPointer();}
        ConstElement_ptr getConstPointer() const { return TransposeMatrix<Matrix>::_Mat.getConstPointer();}

        template <class Vector1, class Vector2>
        Vector1&  apply (Vector1& y, const Vector2& x) const { return TransposeMatrix<Matrix>::_Mat.applyTranpose(y,x);}

        template <class Vector1, class Vector2>
        Vector1&  applyTranpose (Vector1& y, const Vector2& x) const { return TransposeMatrix<Matrix>::_Mat.apply(y,x);}

        void random() { TransposeMatrix<Matrix>::_Mat.random();}
        template<typename RandIter>
        void random(RandIter& G) { TransposeMatrix<Matrix>::_Mat.random(G);}

        Matrix& getMatrix() const { return TransposeMatrix<Matrix>::_Mat;}
	};
    
     template <class T>
    struct isTransposed {
        static const auto value=FFLAS::FflasNoTrans;
    };
    template <class M>
    struct isTransposed<TransposedBlasMatrix<M>> {
        static const auto value=FFLAS::FflasTrans;
    };
    template <class M>
    struct isTransposed<TransposedBlasMatrix<TransposedBlasMatrix<M>>> {
        static const auto value=isTransposed<M>::value;
    };
}


  

#endif // __LINBOX_matrix_densematrix_blas_transposed_matrix_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
