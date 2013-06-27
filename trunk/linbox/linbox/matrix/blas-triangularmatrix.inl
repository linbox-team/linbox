/* linbox/matrix/blas-matrix.h
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *
 * Written by :
 *               Pascal Giorgi  <pascal.giorgi@ens-lyon.fr>
 *               Clément Pernet <clement.pernet@imag.fr>
 *               Brice Boyer    <bboyer@imag.fr>
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

/*!@internal
 * @file matrix/blas-matrix.inl
 * @ingroup matrix
 * A \c BlasMatrix<\c _Field > represents a matrix as an array of
 * <code>_Field</code>s.
 */

#ifndef __LINBOX_blas_triangularmatrix_INL
#define __LINBOX_blas_triangularmatrix_INL

namespace LinBox
{
	template < class _Field, class _Rep >
	TriangularBlasMatrix< _Field, _Rep >::TriangularBlasMatrix (const _Field & F, const size_t m, const size_t n,
							    LINBOX_enum (LinBoxTag::Shape) x,
							    LINBOX_enum (LinBoxTag::Diag) y) :
		BlasMatrix< _Field, _Rep >(F, m, n ) , _uplo(x), _diag(y)
	{}

	template < class _Field, class _Rep >
	TriangularBlasMatrix< _Field, _Rep >::TriangularBlasMatrix (const BlasMatrix< _Field, _Rep >& A,
							    LINBOX_enum (LinBoxTag::Shape) x,
							    LINBOX_enum (LinBoxTag::Diag) y) :
		BlasMatrix< _Field, _Rep >(A) , _uplo(x), _diag(y)
	{}

	template < class _Field, class _Rep >
	TriangularBlasMatrix< _Field, _Rep >::TriangularBlasMatrix (BlasMatrix< _Field, _Rep >& A,
							    LINBOX_enum (LinBoxTag::Shape) x,
							    LINBOX_enum (LinBoxTag::Diag) y) :
		BlasMatrix< _Field, _Rep >(A), _uplo(x), _diag(y)
	{}

	template < class _Field, class _Rep >
	TriangularBlasMatrix< _Field, _Rep >::TriangularBlasMatrix (const TriangularBlasMatrix< _Field, _Rep >& A) :
		BlasMatrix< _Field, _Rep >(A.field(), A.rowdim(),A.coldim()), _uplo(A._uplo), _diag(A._diag)
	{
		switch (A._uplo) {
		case LinBoxTag::Shape::Upper:
			{
				for (size_t i=0;i<A.rowdim();++i)
					for (size_t j=i;j<A.coldim();++j)
						this->setEntry(i,j,A.getEntry(i,j));
				break;
			}
		case LinBoxTag::Shape::Lower:
			{
				for (size_t i=0;i<A.rowdim();++i) {
					for (size_t j=0;j<=i;++j)
						this->setEntry(i,j,A.getEntry(i,j));
				}

				break;
			}
		default:
			throw LinboxError ("Error in copy constructor of TriangularBlasMatrix (incorrect argument)");
		}
	}

	template < class _Field, class _Rep >
	template<class Matrix>
	TriangularBlasMatrix< _Field, _Rep >::TriangularBlasMatrix (const Matrix& A,
							    LINBOX_enum (LinBoxTag::Shape) x,
							    LINBOX_enum (LinBoxTag::Diag) y) :
		BlasMatrix< _Field, _Rep >(A.field(),A.rowdim(),A.coldim()), _uplo(x), _diag(y)
	{
		switch (x) {
		case LinBoxTag::Shape::Upper:
			{
				for (size_t i=0;i<A.rowdim();++i){
					for (size_t j=i;j<A.coldim();++j) {
						Element tmp = A.getEntry(i,j) ;
						this->setEntry(i,j,tmp);
					}
				}
				break;
			}
		case LinBoxTag::Shape::Lower:
			{
				for (size_t i=0;i<A.rowdim();++i) {
					for (size_t j=0;j<=i;++j) {
						Element tmp = A.getEntry(i,j);
						this->setEntry(i,j,tmp);
					}
				}

				break;
			}
		default:
			throw LinboxError ("Error in copy constructor of TriangularBlasMatrix (incorrect argument)");
		}
	}

	template < class _Field, class _Rep >
	LINBOX_enum (LinBoxTag::Shape) TriangularBlasMatrix< _Field, _Rep >::getUpLo() const
	{
		return _uplo;
	}

	template < class _Field, class _Rep >
	LINBOX_enum (LinBoxTag::Diag) TriangularBlasMatrix< _Field, _Rep >::getDiag() const
	{
		return _diag;
	}

}

#endif // __LINBOX_blas_triangularmatrix_INL



// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
