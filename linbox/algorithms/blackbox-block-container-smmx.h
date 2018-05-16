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

#include <algorithm>

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
	template<class _Field1>
	class __BBC_MMHelper {
	private:
		MatrixDomain<_Field1> _MD;
		typedef BlasMatrix<_Field1> Mat;
		
	public:
		__BBC_MMHelper(const _Field1 &F) : _MD(F) {}
		
		inline void mul(Mat &C, const Mat &A, const Mat &B) const {
			_MD.mul(C, A, B);
		}
	};
	
	template<>
	class __BBC_MMHelper<Givaro::Modular<double>> {
	private:
		BlasMatrixDomain<Givaro::Modular<double>> _MD;
		typedef BlasMatrix<Givaro::Modular<double>> Mat;
		
	public:
		__BBC_MMHelper(const Givaro::Modular<double> &F) : _MD(F) {}
		
		inline void mul(Mat &C, const Mat &A, const Mat &B) const {
			_MD.mul(C, A, B);
		}
	};
	
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
		
		Givaro::Timer TW;

	private:
		Field _F;
		__BBC_MMHelper<Field> _MD;
		
		bool left_pre = false;
		FSparseMat _L;
		
		const Blackbox *_BB;
		FSparseMat _M;
		
		bool right_pre = false;
		FSparseMat _R;
		
		Block _U;
		
		size_t _current_data_idx = 0;
		Block _W0;
		Block _W1;
		
		Value _V;
	
		static void convert(const Field &F, FSparseMat &A, const Blackbox *BB, bool transpose = false) {
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
		BlackboxBlockContainerSmmx() {}

		// constructor of the sequence from a blackbox, a field and one block projection
		BlackboxBlockContainerSmmx(
			const Blackbox *BB,
			const Field &F,
			const Block &U0,
			const Block &V0,
			bool transpose = false) : 
				_F(F), 
				_MD(F),
				_BB(BB), 
				_U(U0), 
				_W0(V0), 
				_W1(F, BB->rowdim(), V0.coldim()), 
				_V(F, U0.rowdim(), V0.coldim()) 
		{
			convert(F, _M, BB, transpose);
		}
		
		BlackboxBlockContainerSmmx(
			const Blackbox *BB,
			const Blackbox *PreR,
			const Field &F,
			const Block &U0,
			const Block &V0) : BlackboxBlockContainerSmmx(BB, F, U0, V0)
		{
			convert(F, _R, PreR);
			right_pre = true;
		}
		
		BlackboxBlockContainerSmmx(
			const Blackbox *PreL,
			const Blackbox *BB,
			const Blackbox *PreR,
			const Field &F,
			const Block &U0,
			const Block &V0) : BlackboxBlockContainerSmmx(BB, PreR, F, U0, V0)
		{
			convert(F, _L, PreL);
			left_pre = true;
		}
		
		const Field& field() const {
			return _F;
		}
		
		size_t bbdim() const {
			return _U.coldim();
		}
		
		size_t rowdim() const {
			return _U.rowdim();
		}
		
		size_t coldim() const {
			return _V.coldim();
		}
		
		double spmmTime() const {
			return TW.usertime();
		}
		
		void next() {
			size_t b = coldim();
			
			// W = _R * W
			if (right_pre && _current_data_idx == 0) {
				FFLAS::fspmm(_F, _R, b, _W0.getPointer(), b, _F.zero, _W1.getPointer(), b);
				_current_data_idx = 1;
			} else if (right_pre && _current_data_idx == 1) {
				FFLAS::fspmm(_F, _R, b, _W1.getPointer(), b, _F.zero, _W0.getPointer(), b);
				_current_data_idx = 0;
			}
			
			TW.start();
			// W = _M * W
			if (_current_data_idx == 0) {
				FFLAS::fspmm(_F, _M, b, _W0.getPointer(), b, _F.zero, _W1.getPointer(), b);
				_current_data_idx = 1;
			} else {
				FFLAS::fspmm(_F, _M, b, _W1.getPointer(), b, _F.zero, _W0.getPointer(), b);
				_current_data_idx = 0;
			}
			TW.stop();
				
			// W = _L * W
			if (left_pre && _current_data_idx == 0) {
				FFLAS::fspmm(_F, _L, b, _W0.getPointer(), b, _F.zero, _W1.getPointer(), b);
				_current_data_idx = 1;
			} else if (left_pre && _current_data_idx == 1) {
				FFLAS::fspmm(_F, _L, b, _W1.getPointer(), b, _F.zero, _W0.getPointer(), b);
				_current_data_idx = 0;
			}
		}
		
		const Value &getValue() {
			if (_current_data_idx == 0) {
				_MD.mul(_V, _U, _W0);
			} else {
				_MD.mul(_V, _U, _W1);
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