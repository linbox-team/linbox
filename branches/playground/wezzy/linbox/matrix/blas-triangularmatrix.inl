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
	template<class _Field>
		TriangularBlasMatrix<_Field>::TriangularBlasMatrix (const _Field & F, const size_t m, const size_t n,
				      LinBoxTag::Shape x,
				      LinBoxTag::Diag y) :
			BlasMatrix<_Field>(F, m, n ) , _uplo(x), _diag(y)
		{}

		template<class _Field>
		TriangularBlasMatrix<_Field>::TriangularBlasMatrix (const BlasMatrix<_Field>& A,
				      LinBoxTag::Shape x,
				      LinBoxTag::Diag y) :
			BlasMatrix<_Field>(A) , _uplo(x), _diag(y)
		{}

		template<class _Field>
		TriangularBlasMatrix<_Field>::TriangularBlasMatrix (BlasMatrix<_Field>& A,
				LinBoxTag::Shape x,
				LinBoxTag::Diag y) :
			BlasMatrix<_Field>(A), _uplo(x), _diag(y)
		{}

		template<class _Field>
		TriangularBlasMatrix<_Field>::TriangularBlasMatrix (const TriangularBlasMatrix<_Field>& A) :
			BlasMatrix<_Field>(A.field(), A.rowdim(),A.coldim()), _uplo(A._uplo), _diag(A._diag)
		{
			switch (A._uplo) {
			case LinBoxTag::Upper:
				{
					for (size_t i=0;i<A.rowdim();++i)
						for (size_t j=i;j<A.coldim();++j)
							this->setEntry(i,j,A.getEntry(i,j));
					break;
				}
			case LinBoxTag::Lower:
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

		template<class _Field>
		template<class Matrix>
		TriangularBlasMatrix<_Field>::TriangularBlasMatrix (const Matrix& A,
				LinBoxTag::Shape x,
				LinBoxTag::Diag y) :
			BlasMatrix<_Field>(A.field(),A.rowdim(),A.coldim()), _uplo(x), _diag(y)
		{
			switch (x) {
			case LinBoxTag::Upper:
				{
					for (size_t i=0;i<A.rowdim();++i){
						for (size_t j=i;j<A.coldim();++j) {
							Element tmp = A.getEntry(i,j) ;
							this->setEntry(i,j,tmp);
						}
					}
					break;
				}
			case LinBoxTag::Lower:
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

		template<class _Field>
		LinBoxTag::Shape TriangularBlasMatrix<_Field>::getUpLo() const
		{
			return _uplo;
		}

		template<class _Field>
		LinBoxTag::Diag TriangularBlasMatrix<_Field>::getDiag() const
		{
			return _diag;
		}

}
#endif // __LINBOX_blas_triangularmatrix_INL



// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

