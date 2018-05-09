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

#ifndef __LINBOX_blackbox_block_container_smmx2_H
#define __LINBOX_blackbox_block_container_smmx2_H

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
	class __BBC_MMHelper2 {
	private:
		MatrixDomain<_Field1> _MD;
		typedef BlasMatrix<_Field1> Mat;
		
	public:
		__BBC_MMHelper2(const _Field1 &F) : _MD(F) {}
		
		inline void mul(Mat &C, const Mat &A, const Mat &B) const {
			_MD.mul(C, A, B);
		}
	};
	
	template<>
	class __BBC_MMHelper2<Givaro::Modular<double>> {
	private:
		BlasMatrixDomain<Givaro::Modular<double>> _MD;
		typedef BlasMatrix<Givaro::Modular<double>> Mat;
		
	public:
		__BBC_MMHelper2(const Givaro::Modular<double> &F) : _MD(F) {}
		
		inline void mul(Mat &C, const Mat &A, const Mat &B) const {
			_MD.mul(C, A, B);
		}
	};
	
	template<class _Field>
	class BlackboxBlockContainerSmmx2 {
	public:
		typedef _Field                         Field;
		typedef typename Field::Element        Element;
		typedef typename Field::RandIter       RandIter;
		typedef typename Field::Element_ptr    FflasBlock;
		typedef BlasMatrix<Field>              Block;
		typedef BlasMatrix<Field>              Value;
		      
		typedef FFLAS::Sparse<Field, FFLAS::SparseMatrix_t::CSR> FSparseMat;

	private:
		Field _F;
		__BBC_MMHelper2<Field> _MD;
		
		FSparseMat _M;
		
		Block _U;
		Block _W;
		Block _tmp;
		Value _tmpV;
		Value _V;
		
	public:
		// Default constructor
		BlackboxBlockContainerSmmx2() {}

		// constructor of the sequence from a blackbox, a field and one block projection
		BlackboxBlockContainerSmmx2(
			const FSparseMat &BB,
			const Field &F,
			const BlasVector<Field> &U0,
			const Block &V0) : 
				_F(F), 
				_MD(F),
				_M(BB), 
				_U(F, 1, U0.size()), 
				_W(V0), 
				_tmp(F, U0.size(), V0.coldim()), 
				_tmpV(F, 1, V0.coldim()),
				_V(F, V0.coldim(), 1) 
		{
			for (size_t i = 0; i < U0.size(); i++) {
				_U.setEntry(0, i, U0.getEntry(i));
			}
		}
		
		const Field& field() const {
			return _F;
		}
		
		size_t bbdim() const {
			return _U.coldim();
		}
		
		size_t rowdim() const {
			return _V.rowdim();
		}
		
		size_t coldim() const {
			return _V.coldim();
		}
		
		void next() {
			size_t b = _W.coldim();
			
			// W = _M * W
			for (size_t i = 0; i < _W.rowdim() * _W.coldim(); i++) {
				_tmp.getPointer()[i] = _W.getPointer()[i];
			}
			FFLAS::fspmm(_F, _M, b, _tmp.getPointer(), b, _F.zero, _W.getPointer(), b);
		}
		
		const Value &getValue() {
			_MD.mul(_tmpV, _U, _W);
			for (size_t i = 0; i < _tmpV.coldim(); i++) {
				_V.setEntry(i, 0, _tmpV.getEntry(0, i));
			}
			return _V;
		}
		
		class const_iterator {
			BlackboxBlockContainerSmmx2<Field> *_c;
		public:
			const_iterator() : _c(0){} // BB ??
			const_iterator(BlackboxBlockContainerSmmx2<Field> &C) :_c(&C) {}
			const_iterator &operator++() { _c->next(); return *this; }
			const Value &operator*() { return _c->getValue(); }
		};

		const_iterator begin () { return const_iterator(*this); }
		const_iterator end   () { return const_iterator(); }
	};
}

#endif // __LINBOX_blackbox_block_container_smmx2_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s