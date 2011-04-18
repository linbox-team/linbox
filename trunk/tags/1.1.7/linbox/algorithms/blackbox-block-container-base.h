/* linbox/algorithms/blackbox-block-container-base.h
 * Copyright (C) 2002 Pascal Giorgi
 *
 * Written by Pascal Giorgi pascal.giorgi@ens-lyon.fr
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef __LINBOX_blackbox_block_container_base_H
#define __LINBOX_blackbox_block_container_base_H


#include "time.h" // for seeding

#include <linbox/blackbox/archetype.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/algorithms/blas-domain.h>
#include <linbox/util/debug.h>
#undef _U
#undef _V
#undef _F

namespace LinBox 
{

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

	/** \brief A base class for BlackboxBlockContainer.
	 * The primary member function is begin().

	 * It returns an iterator which after i increments (++) dereferences to 
	 * $U A^i V$, for $U$ and $V$ determined by the init function.
	 * It is designed to be used with implementations of Block Berlekamp-Massey 
	 * such as BlockMasseyDomain.
	 *
	 * Subclasses complete the implementation by defining _launch() and _wait().
	 */

	template<class _Field, class _Blackbox>
	class BlackboxBlockContainerBase {
	public:
		typedef _Field                         Field;
		typedef typename Field::Element      Element;
		typedef _Blackbox                   Blackbox;
		typedef BlasMatrix<Element>            Block;
		typedef BlasMatrix<Element>            Value;

		// Default constructors
		BlackboxBlockContainerBase () {} 

		// Sequence constructor from a blackbox and a field
		// cs set the size of the sequence
		BlackboxBlockContainerBase (const Blackbox *BD, const Field &F, size_t m, size_t n, size_t seed=time(NULL))
			: _F(F)  , _BB(BD), _size(BD->rowdim()/m + BD->coldim()/n +2) , _row(BD->rowdim()),  _col(BD->coldim()), _m(m), _n(n),  _value(m,n), _seed(seed) {}
  
	
		virtual ~BlackboxBlockContainerBase(){}

		// iterator of the sequence
		class const_iterator {
			
		protected:
			BlackboxBlockContainerBase<Field, Blackbox> &_c;

		public:
			const_iterator () {}
			
			const_iterator (BlackboxBlockContainerBase<Field, Blackbox> &C) : _c (C) {}
			
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
		const Field &getField () const { return _F; }

		// blackbox of the sequence
		const Blackbox *getBB () const { return _BB; }		

		// row dimension of the sequence element
		size_t rowdim() const          { return _m; }
		
		// column dimension of the sequence element
		size_t coldim() const          { return _n; }

		// row dimension of the matrix 
		size_t getrow() const { return _BB->rowdim();}

		// col dimension of the matrix 
		size_t getcol() const { return _BB->rowcol();}


	protected:

		friend class const_iterator;
    
		/** Launches a process to do the computation of the next sequence 
		 *  value: $U A^{i+1} V$.  ...or just does it.
		 */
		virtual void _launch() = 0;

		/** If a separate process is computing the next value of $U A^{i+1} V$,
		 * _wait() blocks until the value is ready.
		 */
		virtual void _wait() = 0;

		//-------------- 
		/// Members
		//--------------  
 
		Field                        _F;
		const Blackbox             *_BB;
		size_t                    _size;
		size_t                     _row;
		size_t                     _col;
		size_t                       _m;
		size_t                       _n;
		
                // BDS 22.03.03
		long                 casenumber;
		Block                        _U;
		Block                        _V;
		Value                    _value;
		size_t                    _seed;


		const Value &getvalue() {return _value;}

		//-------------- 
		/// Initializers
		//--------------  

		// Blackbox multiplication using apply function 
		void Mul(Block &M1, const Blackbox &M2 , const Block& M3) {
			linbox_check( M1.rowdim() == M2.rowdim());
			linbox_check( M2.coldim() == M3.rowdim());
			linbox_check( M1.coldim() == M3.coldim());
			
			typename Block::ColIterator        p1 = M1.colBegin();
			typename Block::ConstColIterator   p3 = M3.colBegin();
			
			for (; p3 != M3.colEnd(); ++p1,++p3) {
				M2.apply(*p1,*p3);
			}
		}
	    
	
		/// User Left and Right blocks
		void init (const Block& U, const Block& V) {
			
			linbox_check ( U.coldim() == _row);
			linbox_check ( V.rowdim() == _col);
			casenumber = 1;
			_U = U;  
			_V = V;  
			_value = Value(_m,_n);
			BlasMatrixDomain<Field> BMD(_F);
			BMD.mul(_value, _U, _V);
		}

		// Random Left Matrix and Right Matrix       
		void init (size_t m, size_t n) {
			casenumber = 1;
								
			typename Field::RandIter G(_F,0,_seed);
			Block U (m, _BB->rowdim());
			_U =U;			
			Block V(_BB->coldim(), n);
			_V =V;
			
			typename Block::RawIterator iter_U = _U.rawBegin();
			for (; iter_U != _U.rawEnd();++iter_U)
				G.random(*iter_U);
			
			typename Block::RawIterator iter_V = _V.rawBegin();
			for (; iter_V != _V.rawEnd();++iter_V)
				G.random(*iter_V);

			_value = Value(m,n);
			BlasMatrixDomain<Field> BMD(_F);
			BMD.mul(_value, _U, _V);
		}
	};
 
}

#endif // __LINBOX_blackbox_block_container_base_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
