/* linbox/algorithms/blackbox-block-container.h
 * Copyright (C) 2002 Pascal Giorgi
 *
 * Written by Pascal Giorgi pascal.giorgi@ens-lyon.fr
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

/*! @file algorithms/blackbox-block-container.h
 * @ingroup algorithms
 * @brief no doc.
 */

#ifndef __LINBOX_blackbox_fflas_csr_H
#define __LINBOX_blackbox_fflas_csr_H

#include <algorithm>

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

#include "linbox/blackbox/blockbb.h"

#include "linbox/algorithms/blackbox-block-container-base.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas/fflas_sparse.h"

namespace LinBox
{
	template<class _Field>
	class FflasCsr {
	public:
		typedef _Field                         Field;
		typedef typename Field::Element        Element;
		
		typedef FFLAS::Sparse<Field, FFLAS::SparseMatrix_t::CSR> FSparseMat;
		
	private:
		Field _F;
		FSparseMat _A;
		FSparseMat _T;
		
		size_t _rowdim, _coldim;
	
		static void convert(const Field &F, FSparseMat &A, const SparseMatrix<Field, SparseMatrixFormat::CSR> *BB, bool transpose = false) {
			std::vector<index_t> st1 = BB->getStart();
			std::vector<index_t> col1 = BB->getColid();
			std::vector<Element> data1 = BB->getData();
			
			uint64_t nnz = data1.size();
			
			index_t *row = FFLAS::fflas_new<index_t>(nnz);
			for (size_t j = 0; j < BB->rowdim(); ++j) {
				for (index_t k = st1[j] ; k < st1[j+1]; ++k) {
					row[k] = j;
				}
			}
			
			index_t *col = FFLAS::fflas_new<index_t>(nnz);
			for (size_t i = 0; i < col1.size(); i++) {
				col[i] = col1[i];
			}
			
			typename Field::Element_ptr data = FFLAS::fflas_new<Element>(nnz);
			for (size_t i = 0; i < data1.size(); i++) {
				data[i] = data1[i];
			}
			
			if (transpose) {
				std::swap(row, col);
			}
			
			FFLAS::sparse_init(F, A, row, col, data, BB->rowdim(), BB->coldim(), nnz);
		}
		
	public:
		// Default constructor
		FflasCsr() {}

		// constructor of the sequence from a blackbox, a field and one block projection
		FflasCsr(const SparseMatrix<Field, SparseMatrixFormat::CSR> *BB) : _F(BB->field()) 
		{
			_rowdim = BB->rowdim();
			_coldim = BB->coldim();
			
			convert(_F, _A, BB, false);
			convert(_F, _T, BB, true);
		}
		
		FflasCsr(const FflasCsr &BB) : 
			_F(BB._F), 
			_rowdim(BB._rowdim), 
			_coldim(BB._coldim), 
			_A(BB._A), 
			_T(BB._T) 
		{}
		
		// C = AB
		template<class Matrix>
		Matrix& applyLeft(Matrix &C, const Matrix &B) const {
			size_t b = B.coldim();
			FFLAS::fspmm(_F, _A, b, B.getPointer(), b, _F.zero, C.getPointer(), b);
			return C;
		}
		
		// C = BA
		template<class Matrix>
		Matrix& applyRight(Matrix &C, const Matrix &B) const {
			return C;
		}
	
		template<class OutVector, class InVector>
		OutVector& apply(OutVector& y, const InVector& x) const {
			return y;
		}
		
		template<class OutVector, class InVector>
		OutVector& applyTranspose(OutVector& y, const InVector& x) const {
			return y;
		}
		
		const Field& field() const {
			return _F;
		}
		
		size_t rowdim() const {
			return _rowdim;
		}
		
		size_t coldim() const {
			return _coldim;
		}
	};
	
	template<class Field>
	struct is_blockbb<FflasCsr<Field>> {
		static const bool value = true;
	};
}

#endif // __LINBOX_blackbox_fflas_csr_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
