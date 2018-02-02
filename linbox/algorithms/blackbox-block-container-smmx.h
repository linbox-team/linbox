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

#ifndef __LINBOX_blackbox_block_container_smmx_H
#define __LINBOX_blackbox_block_container_smmx_H

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

#include "linbox/algorithms/blackbox-block-container-base.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas/fflas_sparse.h"

namespace LinBox
{
	template<class _Field, class _Blackbox>
	class BlackboxBlockContainerSmmx {
	public:
		typedef _Field                         Field;
		typedef typename Field::Element        Element;
		typedef typename Field::RandIter       RandIter;
		typedef typename Field::Element_ptr    FflasBlock;
		typedef _Blackbox                      Blackbox;
		typedef BlasMatrix<Field>              Block;
		typedef BlasMatrix<Field>              Value;
		
		typedef FFLAS::Sparse<Field, FFLAS::SparseMatrix_t::CSR> FSparseMat;

	private:
		Field _F;
		const Blackbox *_BB;
		Value _V;
		
		size_t _n;
		size_t _b;
		FSparseMat _M;
		FflasBlock _U;
		FflasBlock _W;
		FflasBlock _C;
		FflasBlock _tmp;
		
	public:
		// Default constructor
		BlackboxBlockContainerSmmx() {}

		// constructor of the sequence from a blackbox, a field and one block projection
		BlackboxBlockContainerSmmx(
			const Blackbox *BB,
			const Field &F,
			const Block &U0,
			const Block &V0) : 
		_F(F), _BB(BB), _V(F, U0.rowdim(), U0.rowdim()) {
			size_t b = U0.rowdim();
			size_t n = BB->rowdim();
			
			_n = n;
			_b = b;
			
			_U = FFLAS::fflas_new(_F, b, n, Alignment::CACHE_PAGESIZE);
			_W = FFLAS::fflas_new(_F, n, b, Alignment::CACHE_LINE);
			_C = FFLAS::fflas_new(_F, b, b, Alignment::CACHE_PAGESIZE);
			_tmp = FFLAS::fflas_new(_F, _n, _b, Alignment::CACHE_LINE);
			
			for (size_t i = 0; i < b; i++) {
				for (size_t j = 0; j < n; j++) {
					_U[i * n + j] = U0.getEntry(i, j);
				}
			}
			
			for (size_t i = 0; i < n; i++) {
				for (size_t j = 0; j < b; j++) {
					_W[i * b + j] = V0.getEntry(i, j);
				}
			}
			
			Blackbox M(*BB);
			std::vector<index_t> st1 = M.getStart();
			std::vector<index_t> col1 = M.getColid();
			std::vector<Element> data1 = M.getData();
			
			uint64_t nnz = data1.size();
			
			index_t *row = FFLAS::fflas_new<index_t>(nnz);
			for (size_t j = 0; j < M.rowdim(); ++j) {
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
			
			FFLAS::sparse_init(F, _M, row, col, data, M.rowdim(), M.coldim(), nnz);
		}
		
		const Field& field() const {
			return _F;
		}
		
		const Blackbox *getBB() const {
			return _BB;
		}
		
		size_t rowdim() const {
			return _b;
		}
		
		size_t coldim() const {
			return _b;
		}
		
		void next() {
			for (size_t i = 0; i < _n * _b; i++) _tmp[i] = _W[i];
			FFLAS::fspmm(_F, _M, _b, _tmp, _b, _F.zero, _W, _b);
		}
		
		const Value &getValue() {
			fgemm(_F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, _b, _b, _n, _F.one, _U, _n, _W, _b, _F.zero, _C, _b);
			
			for (size_t i = 0; i < _b; i++) {
				for (size_t j = 0; j < _b; j++) {
					_V.setEntry(i, j, _C[i * _b + j]);
				}
			}
			
			return _V;
		}
		
		class const_iterator {
			BlackboxBlockContainerSmmx<Field, Blackbox> *_c;
		public:
			const_iterator() : _c(0){} // BB ??
			const_iterator(BlackboxBlockContainerSmmx<Field, Blackbox> &C) :_c(&C) {}
			const_iterator &operator++() { _c->next(); return *this; }
			const Value &operator*() { return _c->getValue(); }
		};

		const_iterator begin () { return const_iterator(*this); }
		const_iterator end   () { return const_iterator(); }
	};
}

#endif // __LINBOX_blackbox_block_container_smmx_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s