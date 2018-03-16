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

#ifndef __LINBOX_blackbox_block_container_spmv_H
#define __LINBOX_blackbox_block_container_spmv_H

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

#include "linbox/algorithms/blackbox-block-container-base.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "fflas-ffpack/config-blas.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/fflas/fflas_sparse.h"

#include "givaro/givtimer.h"

namespace LinBox
{
	template<class _Field, class _Blackbox>
	class BlackboxBlockContainerSpmv {
	public:
		typedef _Field                         Field;
		typedef typename Field::Element        Element;
		typedef typename Field::RandIter       RandIter;
		typedef typename Field::Element_ptr    FflasBlock;
		typedef _Blackbox                      Blackbox;
		typedef BlasMatrix<Field>              Block;
		typedef BlasMatrix<Field>              Value;
		
	private:
		Field _F;
		BlasMatrixDomain<Field> _BMD;
		const Blackbox *_BB;
		Block _U;
		Block _W;
		Block _tmp;
		Block _V;
		
		Givaro::Timer TW;
		double _spmv_time;
		double _gemm_time;
		
	public:
		// Default constructor
		BlackboxBlockContainerSpmv() {}

		// constructor of the sequence from a blackbox, a field and one block projection
		BlackboxBlockContainerSpmv(
			const Blackbox *BB,
			const Field &F,
			const Block &U0,
			const Block &V0) : 
		_F(F), _BMD(F), _BB(BB), _U(U0), _W(V0), _tmp(V0), _V(F, U0.rowdim(), U0.rowdim()) {
			_spmv_time = 0;
			_gemm_time = 0;
		}
		
		const Field& field() const {
			return _F;
		}
		
		const Blackbox *getBB() const {
			return _BB;
		}
		
		size_t rowdim() const {
			return _V.rowdim();
		}
		
		size_t coldim() const {
			return _V.coldim();
		}
		
		double spmmtime() const {
			return _spmv_time;
		}
		
		double gemmtime() const {
			return _gemm_time;
		}
		
		void next() {
			MatrixDomain<Field> MD(_F);
			typename Block::ColIterator p1 = _tmp.colBegin();
			typename Block::ConstColIterator p2 = _W.colBegin();
			
			TW.clear();
			TW.start();
			
			for (; p2 != _W.colEnd(); ++p1, ++p2) {
				_BB->apply(*p1, *p2);
			}
			
			MD.copy(_W, _tmp);
			
			TW.stop();
			_spmv_time += TW.usertime();
		}
		
		const Value &getValue() {
			TW.clear();
			TW.start();
			
			_BMD.mul(_V, _U, _W);
			
			TW.stop();
			_gemm_time += TW.usertime();
			
			return _V;
		}
		
		class const_iterator {
			BlackboxBlockContainerSpmv<Field, Blackbox> *_c;
		public:
			const_iterator() : _c(0){} // BB ??
			const_iterator(BlackboxBlockContainerSpmv<Field, Blackbox> &C) :_c(&C) {}
			const_iterator &operator++() { _c->next(); return *this; }
			const Value &operator*() { return _c->getValue(); }
		};

		const_iterator begin () { return const_iterator(*this); }
		const_iterator end   () { return const_iterator(); }
	};
}

#endif // __LINBOX_blackbox_block_container_spmv_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s