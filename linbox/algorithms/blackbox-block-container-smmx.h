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
		typedef typename Field::Element_ptr    Block;
		typedef _Blackbox                      Blackbox;
		typedef BlasMatrix<Field>              Value;

	private:
		Field _F;
		const Blackbox *_M;
		
		Block _U;
		Block _W;
		size_t _b;
		Value _V;
		
		Block _C;
		
	public:
		// Default constructor
		BlackboxBlockContainerSmmx() {}

		// constructor of the sequence from a blackbox, a field and one block projection
		BlackboxBlockContainerSmmx(
			const _Blackbox *M,
			const Field &F,
			const Block &U0,
			const Block& V0,
			size_t b) : 
		_F(F), _M(M), _U(U0), _W(V0), _b(b), _V(_F, b, b) {
			_C = FFLAS::fflas_new(_F, b, b, Alignment::CACHE_PAGESIZE);
			
			int nbw = -1;
			FFLAS::MMHelper<Field, FFLAS::MMHelperAlgo::Winograd> WH (F, nbw, FFLAS::ParSeqHelper::Sequential());
			fgemm(_F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, b, b, _M->m, F.one, _U, _M->m, _W, b, F.zero, _C, b, WH);
		}
		
		const Field& field() const {
			return _F;
		}
		
		const _Blackbox *getBB() const {
			return _M;
		}
		
		size_t bbrowdim() const {
			return _M->m;
		}
		
		size_t rowdim() const {
			return _b;
		}
		
		size_t coldim() const {
			return _b;
		}
		
		void next() {
			Block tmp(_W);
			
			FFLAS::fspmm(_F, *_M, coldim(), tmp, coldim(), 1, _W, coldim());
			
			int nbw = -1;
			FFLAS::MMHelper<Field, FFLAS::MMHelperAlgo::Winograd> WH (_F, nbw, FFLAS::ParSeqHelper::Sequential());
			fgemm(_F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, _b, _b, _M->m, _F.one, _U, _M->m, _W, _b, _F.zero, _C, _b, WH);
			
			for (size_t i = 0; i < _b * _b; i++) {
				_V.setEntry(i / _b, i % _b, _C[i]);
			}
		}
		
		const Value &getValue() const {
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
