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

#ifndef __LINBOX_blackbox_block_container_H
#define __LINBOX_blackbox_block_container_H

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

#include "linbox/algorithms/blackbox-block-container-base.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#define _BBC_TIMING

#ifdef _BBC_TIMING
#include <time.h>
#include "linbox/util/timer.h"
#endif

namespace LinBox
{

	// For non-block blackbox's use blockbb to convert to a block blackbox
	template<class _Field, class _Blackbox, class _MatrixDomain = BlasMatrixDomain<_Field>>
	class BlackboxBlockContainer : public BlackboxBlockContainerBase<_Field,_Blackbox,_MatrixDomain> {
	public:
		typedef _Field                         Field;
		typedef typename Field::Element      Element;
		typedef typename Field::RandIter   RandIter;
		typedef BlasMatrix<Field>           Block;
		typedef BlasMatrix<Field>           Value;

		// Default constructor
		BlackboxBlockContainer () {}

		// constructor of the sequence from a blackbox, a field and one block projection
		BlackboxBlockContainer(const _Blackbox *D, const Field &F, const Block  &U0) :
			BlackboxBlockContainerBase<Field,_Blackbox,_MatrixDomain> (D, F, U0.rowdim(), U0.coldim()) , _blockW(D->rowdim(), U0.coldim()), _BMD(F)
		{
#ifdef _BBC_TIMING
			clearTimer();
			tSequence.clear();
			tSequence.start();
#endif
			this->init (U0, U0);
#ifdef _BBC_TIMING
			tSequence.stop();
			ttSequence += tSequence;
#endif
		}

		// constructor of the sequence from a blackbox, a field and two blocks projection
		BlackboxBlockContainer(const _Blackbox *D, const Field &F, const Block &U0, const Block& V0) :
			BlackboxBlockContainerBase<Field,_Blackbox,_MatrixDomain> (D, F,U0.rowdim(), V0.coldim())
			, _blockW(F,D->rowdim(), V0.coldim()), _BMD(F)
		{
#ifdef _BBC_TIMING
			clearTimer();
			tSequence.clear();
			tSequence.start();
#endif
			this->init (U0, V0);
#ifdef _BBC_TIMING
			tSequence.stop();
			ttSequence += tSequence;
#endif
		}

		//  constructor of the sequence from a blackbox, a field and two blocks random projection
            BlackboxBlockContainer(const _Blackbox *D, const Field &F, size_t m, size_t n, size_t seed= static_cast<size_t>(std::time(nullptr))) :
			BlackboxBlockContainerBase<Field, _Blackbox, _MatrixDomain> (D, F, m, n,seed)
			, _blockW(F,D->rowdim(), n), _BMD(F)
		{
#ifdef _BBC_TIMING
			clearTimer();
			tSequence.clear();
			tSequence.start();
#endif
			this->init (m, n);
#ifdef _BBC_TIMING
			tSequence.stop();
			ttSequence += tSequence;
#endif
		}


#ifdef _BBC_TIMING
		void clearTimer() {
			ttSequence.clear();
		}

		void printTimer() {
			std::cout<<"Sequence Computation "<<ttSequence<<std::endl<<std::endl;
		}
#endif

	protected:
		Block                        _blockW;
		_MatrixDomain    _BMD;

#ifdef _BBC_TIMING
		Timer     ttSequence, tSequence;
#endif


		// launcher of the next sequence element computation
		void _launch () {
#ifdef _BBC_TIMING
			tSequence.clear();
			tSequence.start();
#endif
			if (this->casenumber) {
                                this->Mul(_blockW,*this->_BB,this->_blockV);
				_BMD.mul(this->_value, this->_blockU, _blockW);
				this->casenumber = 0;
                        }
			else {
                                this->Mul(this->_blockV,*this->_BB,_blockW);
				_BMD.mul(this->_value, this->_blockU, this->_blockV);
				this->casenumber = 1;
			}
#ifdef _BBC_TIMING
			tSequence.stop();
			ttSequence +=tSequence;
#endif
		}

		void _wait () {}
	};

	/*! @brief no doc.
	 */
	template<class _Field, class _Blackbox, class _MatrixDomain = MatrixDomain<_Field>>
	class BlackboxBlockContainerRecord : public BlackboxBlockContainerBase<_Field,_Blackbox,_MatrixDomain> {

	public:
		typedef _Field                        Field;
		typedef typename Field::Element     Element;
		typedef typename Field::RandIter   RandIter;
		typedef BlasMatrix<Field>           Block;
		typedef BlasMatrix<Field>           Value;

		enum Launcher {RowUpdate=0, ColUpdate=1, Nothing=2};

		// Default constructor
		BlackboxBlockContainerRecord () {}

		// constructor of the sequence from a blackbox, a field and one block projection
		BlackboxBlockContainerRecord(const _Blackbox *D, const Field &F, const Block  &U0) :
			BlackboxBlockContainerBase<Field,_Blackbox,_MatrixDomain> (D, F, U0.rowdim(), U0.coldim())
			, _blockW(F, D->rowdim(), U0.coldim()), _BMD(F), _launcher(Nothing), _iter(1)
		{
			this->init (U0, U0);
			_rep = std::vector<Value> (this->_size, Value(F));
			_Vcopy = this->_blockV;
			for (size_t i=0;i< this->_size;++i){
				_rep[i] = this->_value;
				_launch_record();
			}
			this->_value=_rep[0];
		}

		// constructor of the sequence from a blackbox, a field and two blocks projection
		BlackboxBlockContainerRecord(const _Blackbox *D, const Field &F, const Block &U0, const Block& V0, bool denseblock= true) :
			BlackboxBlockContainerBase<Field,_Blackbox,_MatrixDomain> (D, F,U0.rowdim(), V0.coldim())
			, _blockW(F,D->rowdim(), V0.coldim()), _BMD(F),  _launcher(Nothing), _iter(1)
		{
#ifdef _BBC_TIMING
			clearTimer();
			tSequence.clear();
			tSequence.start();
#endif
			this->init (U0, V0);


			_rep = std::vector<Value> (this->_size, Value(F));
			_Vcopy = this->_blockV;
			if (denseblock) {
				for (size_t i=0;i< this->_size;++i){
					_rep[i] = this->_value;
					_launch_record();
				}
			}
			else {
				_rep.resize(_rep.size()-3);
				size_t block= U0.rowdim();
				size_t Nn = U0.coldim()/block;
				_Special_U.resize(block, std::vector<typename Field::Element>(Nn));

				for (size_t i=0;i<block;++i)
					for (size_t j=0;j<Nn;++j)
						F.assign(_Special_U[i][j], U0.getEntry(0, i*Nn+j));

				for (size_t i=0;i< this->_size-3;++i){
					_launch_record_notdense();
					_rep[i] = this->_value;
				}

			}

			this->_value=_rep[0];
#ifdef _BBC_TIMING
			tSequence.stop();
			ttSequence += tSequence;
#endif
		}

		//  constructor of the sequence from a blackbox, a field and two blocks random projection
            BlackboxBlockContainerRecord(const _Blackbox *D, const Field &F, size_t m, size_t n, size_t seed= static_cast<size_t>(std::time(nullptr))) :
			BlackboxBlockContainerBase<Field, _Blackbox, _MatrixDomain> (D, F, m, n,seed),
			_blockW(D->rowdim(), n), _BMD(F), _launcher(Nothing), _iter(1)
		{
#ifdef _BBC_TIMING
			clearTimer();
			tSequence.clear();
			tSequence.start();
#endif
			this->init (m,n);
			_rep = std::vector<Value> (this->_size);
			_Vcopy = this->_blockV;
			for (size_t i=0;i< this->_size;++i){
				_rep[i] = this->_value;
				_launch_record();
			}
			this->_value=_rep[0];
#ifdef _BBC_TIMING
			tSequence.stop();
			ttSequence += tSequence;
#endif
		}


		void setU (const std::vector<Element> &b, size_t k)
		{
			linbox_check( b.size() == this->_row);
			_launcher = RowUpdate;
			_upd_idx  = k;
			_u        = b;
			_w.resize(this->_row);
			_iter     = 1;
			_case     = 1;
#ifdef _BBC_TIMING
			tSequence.clear();
			tSequence.start();
#endif
			std::vector<Element> _row_value(this->_n);
			_BMD.mul(_row_value, b, _Vcopy);
			this->_value  = _rep[0];
			for (size_t j=0; j< this->_n; ++j)
				this->_value.setEntry(_upd_idx, j, _row_value[j]);
#ifdef _BBC_TIMING
			tSequence.stop();
			ttSequence += tSequence;
#endif
		}

		void setV (const std::vector<Element> &b, size_t k)
		{

			linbox_check( b.size() == this->_col);
			_launcher = ColUpdate;
			_upd_idx  = k;
			_u        = b;
			_w.resize(this->_col);
			_iter     = 1;
			_case     = 1;
#ifdef _BBC_TIMING
			tSequence.clear();
			tSequence.start();
#endif
			std::vector<Element> _col_value(this->_m);
			_BMD.mul(_col_value, this->_blockU, b);
			this->_value  = _rep[0];
			for (size_t j=0; j< this->_m; ++j)
				this->_value.setEntry(j, _upd_idx, _col_value[j]);
#ifdef __BBC_TIMING
			tSequence.stop();
			ttSequence += tSequence;
#endif
		}


		void recompute()
		{
			switch(_launcher) {
			case Nothing:
				this->_value=_rep[0];
				_iter=1;
				break;
			case RowUpdate:
#ifdef _BBC_TIMING
				tSequence.clear();
				tSequence.start();
#endif
				for (size_t i=0;i< this->_size;++i){
					_rep[i]=this->_value;
					_launch_record_row();
				}
				_launcher=Nothing;
				this->_value=_rep[0];
				_iter=1;
#ifdef _BBC_TIMING
				tSequence.stop();
				ttSequence += tSequence;
#endif
				break;
			case ColUpdate:
#ifdef _BBC_TIMING
				tSequence.clear();
				tSequence.start();
#endif
				for (size_t i=0;i< this->_size;++i){
					_rep[i]=this->_value;
					_launch_record_col();
				}
				_launcher=Nothing;
				this->_value=_rep[0];
				_iter=1;
				break;
#ifdef _BBC_TIMING
				tSequence.stop();
				ttSequence += tSequence;
#endif
			default :
				throw LinboxError ("Bad argument in BlackboxBlockContainerRecord, _launch() function\n");
				break;
			}
		}


#ifdef _BBC_TIMING
		void clearTimer()
		{
			ttSequence.clear();
		}

		void printTimer()
		{
			std::cout<<"Sequence Computation "<<ttSequence<<std::endl<<std::endl;
		}
#endif

		const std::vector<Value>& getRep() const { return _rep;}


	protected:

		Block                           _blockW;
		Block                       _Vcopy;
		_MatrixDomain       _BMD;
		std::vector<Value>            _rep;
		size_t                    _upd_idx;
		std::vector<Element>            _u;
		std::vector<Element>            _w;
		Launcher                 _launcher;
		size_t                       _iter;
		size_t                       _case;
		std::vector<std::vector<Element> > _Special_U;
#ifdef _BBC_TIMING
		Timer        ttSequence, tSequence;
#endif

		// launcher of computation of sequence element
		void _launch_record ()
		{
			if (this->casenumber) {
				Mul(_blockW,*this->_BB,this->_blockV);
				_BMD.mul(this->_value, this->_blockU, _blockW);
				this->casenumber = 0;
			}
			else {
				Mul(this->_blockV,*this->_BB,_blockW);
				_BMD.mul(this->_value, this->_blockU, this->_blockV);
				this->casenumber = 1;
			}
		}


		// launcher of computation of sequence element
		void _launch_record_notdense ()
		{
			size_t block= this->_blockV.coldim();
			size_t numblock=_Special_U[0].size();
			if (this->casenumber) {
				Mul(_blockW,*this->_BB,this->_blockV);

				std::vector<Element> tmp(block);
				for (size_t i=0; i<block; ++i){
					BlasMatrix<Field> T(_blockW, i*numblock, 0, numblock, block);
					_BMD.mul(tmp, _Special_U[i], T);
					for (size_t j=0;j<block;++j){
						this->getField()->assign(this->_value.refEntry(i,j), tmp[j]);
					}
				}

				this->casenumber = 0;
			}
			else {
				Mul(this->_blockV,*this->_BB,_blockW);

				std::vector<Element> tmp(block);
				for (size_t i=0; i< block; ++i){
					BlasMatrix<Field> T(this->_blockV, i*numblock, 0, numblock, block);
					_BMD.mul(tmp, _Special_U[i], T);
					for (size_t j=0;j<block;++j)
						this->getField().assign(this->_value.refEntry(i,j), tmp[j]);
				}

				this->casenumber = 1;
			}
		}


		// launcher of computation of sequence element
		// just update one row in the sequence element
		void _launch_record_row()
		{
			if ( _iter < this->_size) {
				if ( _case == 1) {
					this->_BB->applyTranspose(_w,_u);
					std::vector<Element> _row_value(this->_n);
					_BMD.mul(_row_value, _w, _Vcopy);
					this->_value  = _rep[_iter];
					for (size_t j=0; j< this->_n; ++j)
						this->_value.setEntry(_upd_idx, j, _row_value[j]);
					++_iter;
					_case =0;
				}
				else {
					this->_BB->applyTranspose(_u,_w);
					std::vector<Element> _row_value(this->_n);
					_BMD.mul(_row_value, _u, _Vcopy);
					this->_value  = _rep[_iter];
					for (size_t j=0; j< this->_n; ++j)
						this->_value.setEntry(_upd_idx, j, _row_value[j]);
					++_iter;
					_case =1;
				}
			}
		}

		// launcher of computation of sequence element
		// just update one col in the sequence element
		void _launch_record_col()
		{
			if ( _iter < this->_size) {
				if ( _case == 1) {
					this->_BB->apply(_w,_u);
					std::vector<Element> _col_value(this->_m);
					_BMD.mul(_col_value, this->_blockU, _w);
					this->_value  = _rep[_iter];
					for (size_t j=0; j< this->_m; ++j)
						this->_value.setEntry(j, _upd_idx, _col_value[j]);
					++_iter;
					_case =0;
				}
				else {
					this->_BB->apply(_u,_w);
					std::vector<Element> _col_value(this->_m);
					_BMD.mul(_col_value, this->_blockU, _u);
					this->_value  = _rep[_iter];
					for (size_t j=0; j< this->_m; ++j)
						this->_value.setEntry(j, _upd_idx, _col_value[j]);
					++_iter;
					_case =1;
				}
			}
		}

		// launcher which be used as iterator on the sequence
		void _launch()
		{
			switch (_launcher) {
			case Nothing:
				if (_iter < this->_size)
					this->_value = _rep[_iter];
				++_iter;
				break;
			case RowUpdate:
				_launch_record_row();
				break;
			case ColUpdate:
				_launch_record_col();
				break;
			default :
				throw LinboxError ("Bad argument in BlackboxBlockContainerRecord, _launch() function\n");
				break;
			}
		}

		void _wait () {}
	};

}

#undef _BBC_TIMING

#endif // __LINBOX_blackbox_block_container_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
