/* linbox/algorithms/blackbox-block-container-base.h
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

/*! @file algorithms/blackbox-block-container-base.h
 * @ingroup algorithms
 * @brief NO DOC
 */



#ifndef __LINBOX_blackbox_block_container_base_H
#define __LINBOX_blackbox_block_container_base_H


#include <time.h> // for seeding

#include <omp.h>

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/blackbox/archetype.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/matrix/matrix-domain.h"

// #include "linbox/blackbox/triplesbb-omp.h"
#include "linbox/matrix/sparse-matrix.h"

namespace LinBox
{

//Temporary fix to deal with the fact that not all Blackboxes have applyLeft()
template<class Field,class Block>
class MulHelper {
public:
	template<class Blackbox>
	static void mul(const Field& F,
			Block &M1, const Blackbox &M2, const Block& M3) {
		linbox_check( M1.rowdim() == M2.rowdim());
		linbox_check( M2.coldim() == M3.rowdim());
		linbox_check( M1.coldim() == M3.coldim());

		MatrixDomain<Field> MD(F);
		typename Block::ColIterator        p1 = M1.colBegin();
		typename Block::ConstColIterator   p3 = M3.colBegin();

		for (; p3 != M3.colEnd(); ++p1,++p3) {
			M2.apply(*p1,*p3);
		}
	}

	static void mul (const Field& F,
			 Block &M1, const SparseMatrix<Field,SparseMatrixFormat::TPL> &M2, const Block& M3) {
		M2.applyLeft(M1,M3);
	}
};

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

	/** \brief A base class for BlackboxBlockContainer.
	 * The primary member function is \c begin().

	 * It returns an iterator which after i increments (++) dereferences to
	 * \f$U A^i V\f$, for \f$U\f$ and \f$V\f$ determined by the init
	 * function.  It is designed to be used with implementations of Block
	 * Berlekamp-Massey such as BlockMasseyDomain.
	 *
	 * Subclasses complete the implementation by defining \c _launch() and
	 * \c _wait().
	 */
	template<class _Field, class _Blackbox>
	class BlackboxBlockContainerBase {
	public:
		typedef _Field                         Field;
		typedef typename Field::Element      Element;
		typedef _Blackbox                   Blackbox;
		typedef BlasMatrix<Field>              Block;
		typedef BlasMatrix<Field>              Value;

		// Default constructors
		BlackboxBlockContainerBase () {}

		// Sequence constructor from a blackbox and a field
		// cs set the size of the sequence
		BlackboxBlockContainerBase (const Blackbox *BD, const Field &F, size_t m, size_t n, size_t seed=(size_t)time(NULL)) :
			_field(&F)  , _BB(BD), _size(BD->rowdim()/m + BD->coldim()/n +2)
			, _nn(BD->rowdim()),  _m(m), _n(n)
			,casenumber(0)
			,_blockU(F,_m,_nn),_blockV(F,_nn,_n),_value(field(),m,n), _seed(seed)
		{}


		virtual ~BlackboxBlockContainerBase(){}

		// iterator of the sequence
		class const_iterator {

		protected:
			BlackboxBlockContainerBase<Field, Blackbox> &_c;

		public:
			const_iterator () {}

			const_iterator (BlackboxBlockContainerBase<Field, Blackbox> &C) :
				_c (C)
			{}

			const_iterator &operator ++ () { _c._launch (); return *this; }

			const Value    &operator * ()  { _c._wait (); return _c.getvalue();}
		};

		// begin of the sequence iterator
		const_iterator begin ()        { return const_iterator (*this); }

		// end of the sequence iterator
		const_iterator end ()          { return const_iterator (); }

		// size of the sequence
		size_t size() const            { return _size; }

		// field of the sequence
		const Field &field () const { return *_field; }
		const Field &getField () const { return *_field; } // deprecated

		// blackbox of the sequence
		const Blackbox *getBB () const { return _BB; }

		// row dimension of the sequence element
		size_t rowdim() const          { return _m; }

		// column dimension of the sequence element
		size_t coldim() const          { return _n; }

	protected:

		friend class const_iterator;

		/** Launches a process to do the computation of the next sequence
		 *  value: \f$U A^{i+1} V\f$.  ...or just does it.
		 */
		virtual void _launch() = 0;

		/** If a separate process is computing the next value of \f$U
		 * A^{i+1} V\f$, \c _wait() blocks until the value is ready.
		 */
		virtual void _wait() = 0;

		//--------------
		/// Members
		//--------------

		const Field                        *_field;
		const Blackbox             *_BB;
		size_t                    _size; // length of sequence
		size_t			     _nn; // _BB order (square mat)
		size_t                       _m; // block rows
		size_t                       _n; // block cols

		// BDS 22.03.03
		long                 casenumber;
		Block                        _blockU;
		Block                        _blockV;
		Value                    _value;
		size_t                    _seed;


		const Value &getvalue() {return _value;}

		//--------------
		/// Initializers
		//--------------


		// Blackbox multiplication using apply function
		inline void Mul(Block &M1, const Blackbox &M2 , const Block& M3)
		{
                        MulHelper<Field,Block>::mul(field(),M1,M2,M3);
                }

		/// User Left and Right blocks
		void init (const Block& U, const Block& V)
		{

			linbox_check ( U.rowdim() == _m);
			linbox_check ( U.coldim() == _nn);
			linbox_check ( V.rowdim() == _nn);
			linbox_check ( V.coldim() == _n);
			casenumber = 1;
			_blockU = U;
			_blockV = V;
			_value = Value(*_field,_m,_n);
			BlasMatrixDomain<Field> BMD(*_field);
			BMD.mul(_value, _blockU, _blockV);
		}

		// Random Left Matrix and Right Matrix
		void init (size_t m, size_t n)
		{
			casenumber = 1;

			typename Field::RandIter G(*_field,0,_seed);
			Block U (m, _BB->rowdim());
			_blockU =U;
			Block V(_BB->coldim(), n);
			_blockV =V;

			typename Block::Iterator iter_U = _blockU.Begin();
			for (; iter_U != _blockU.End();++iter_U)
				G.random(*iter_U);

			typename Block::Iterator iter_V = _blockV.Begin();
			for (; iter_V != _blockV.End();++iter_V)
				G.random(*iter_V);

			_value = Value(*_field,m,n);
			BlasMatrixDomain<Field> BMD(*_field);
			BMD.mul(_value, _blockU, _blockV);
		}
	};

}

#endif // __LINBOX_blackbox_block_container_base_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
