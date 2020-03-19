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


namespace LinBox
{ /*  Transposed Matrix */
	/*! TransposedBlasMatrix.
	 * NO DOC
	 */
	template< class Matrix >
	class TransposedBlasMatrix {
	public:
        typedef TransposedBlasMatrix<Matrix> Self_t;
        typedef typename Matrix::Field                    Field;
        typedef typename Field::Element                 Element;  
        typedef typename Field::Element_ptr         Element_ptr;      //!< Pointer to Element type
        typedef typename Field::ConstElement_ptr   ConstElement_ptr; //!< Pointer to const Element type
        

		TransposedBlasMatrix ( Matrix& Mat ) :	_Mat(Mat) {}

        size_t rowdim() const { return _Mat.coldim();}

        size_t coldim() const { return _Mat.rowdim();}

        Self_t& copy(const Self_t& A) { return *this=A;}

        size_t getStride() const { return _Mat.getStride();}
        
        ConstElement_ptr getPointer() const { return _Mat.getPointer();}
        Element_ptr getPointer() { return _Mat.getPointer();}
        ConstElement_ptr getConstPointer() const { return _Mat.getConstPointer();}

        void setEntry (size_t i, size_t j, const Element &a_ij) { _Mat.setEntry(j,i,a_ij);}
        Element&  refEntry (size_t i, size_t j) {return _Mat.refEntry(j,i);}
        Element&  getfEntry (size_t i, size_t j) const {return _Mat.getEntry(j,i);}
        Element&  getfEntry (Element& x, size_t i, size_t j) const {return _Mat.getEntry(x,j,i);}
        const Field& field() const {return _Mat.field();}

        template <class Vector1, class Vector2>
        Vector1&  apply (Vector1& y, const Vector2& x) const { return _Mat.applyTranpose(y,x);}

        template <class Vector1, class Vector2>
        Vector1&  applyTranpose (Vector1& y, const Vector2& x) const { return _Mat.apply(y,x);}

        void random() { _Mat.random();}
        template<typename RandIter>
        void random(RandIter& G) { _Mat.random(G);}

        Matrix& getMatrix() const { return _Mat;}
	protected:
		Matrix& _Mat; //!< NO DOC
	};

	/*! TransposedBlasMatrix.
	 * NO DOC
	 */
	template< class Matrix >
	class TransposedBlasMatrix< TransposedBlasMatrix< Matrix > > : public Matrix {

	public:
		/*! TransposedBlasMatrix.
		 * NO DOC
		 */
		TransposedBlasMatrix ( Matrix& Mat ) :
			Matrix(Mat)
		{}

		/*! TransposedBlasMatrix.
		 * NO DOC
		 */
		TransposedBlasMatrix ( const Matrix& Mat ) :
			Matrix(Mat)
		{}
	protected:
		Matrix& _Mat; //!< NO DOC
        
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
