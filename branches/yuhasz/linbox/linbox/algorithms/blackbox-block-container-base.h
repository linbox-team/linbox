/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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


#ifndef __BLACKBOX_BLOCK_CONTAINER_BASE_H
#define __BLACKBOX_BLOCK_CONTAINER_BASE_H


#include "linbox/blackbox/archetype.h"
#include "linbox/vector/vector-domain.h"
namespace LinBox 
{

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

	/** @memo A base class for BlackboxBlockContainer.
	 * The primary member function is begin().
	 * @doc
	 * It returns an iterator which after i increments (++) dereferences to 
	 * $u A^i v$, for $u$ and $v$ determined by the form of construction.
	 * It is designed to be used with implementations of Berlekamp-Massey by block
	 * such as MasseyBlockDomain.
	 *
	 * Subclasses complete the implementation by defining _launch() and _wait().
	 */

	template<class Field, class Vector,class Matrix>
	class BlackboxBlockContainerBase {
	public:
		typedef BlackboxArchetype<Vector > Blackbox;
		typedef typename Field::Element Element;
		// this Block type should be templatized but I don't know 
		// how to do it in the best way.
		typedef std::vector<Vector> Block;

		//-- Constructors
		BlackboxBlockContainerBase () {} 

		BlackboxBlockContainerBase (const Blackbox *BD, const Field &F)
			: _field (F) , _VD(F) , _BB(BD->clone()), _size (MIN (BD->rowdim (), BD->coldim ()) << 1) , _value(F) { _field.init(Zero,0UL);}
  
	
  
		class const_iterator {
			BlackboxBlockContainerBase<Field, Vector,Matrix> &_c;
		public:
			const_iterator () {}
			const_iterator (BlackboxBlockContainerBase<Field, Vector,Matrix> &C) : _c (C) {}
			const_iterator &operator ++ () { _c._launch (); return *this; }
	
			const Matrix &operator * () { _c._wait (); return _c.getvalue();}			};

		const_iterator begin () { return const_iterator (*this); }
		const_iterator end () { return const_iterator (); }

		long         size     () const { return _size; }
		const Field &getField () const { return _field; }
		Blackbox *getBB () const { return _BB->clone(); }
		size_t row() const {return _row;}
		size_t col() const {return _col;}


	protected:

		friend class const_iterator;
    
		/** Launches a process to do the computation of the next sequence 
		 *  value: $v^T A^{i+1} u$.  ...or just does it.
		 */
		virtual void _launch() = 0;

		/** If a separate process is computing the next value of $v^T A^{i+1} u$,
		 * _wait() blocks until the value is ready.
		 */
		virtual void _wait() = 0;

		//-------------- 
		/// Members
		//--------------  
 
		Field _field;
		Element Zero;
		Blackbox *_BB;
		VectorDomain<Field> _VD;
		long _size;

		size_t _row;
		size_t _col;

                // BDS 22.03.03
		long casenumber;
		Block u;
		Block v;
		Matrix _value;

		const Matrix &getvalue() {return _value;}

		//-------------- 
		/// Initializers
		//--------------  

		void Mul(Matrix &M1, const Block& M2 , const Block& M3) {
			Element Elt;
			_field.init(Elt,0UL);
			long i=0;
			typename Block::const_iterator p2,p3;

			for (p2 = M2.begin(); p2 != M2.end();p2++) {
				int j=0;
				for (p3 = M3.begin(); p3 != M3.end();p3++) {
					_VD.dot(Elt,*p2,*p3);
					M1.setEntry(i,j,Elt);
					j++;
				}
				i++;
			}
		
		}
		void Mul(Block &M1, const Blackbox &M2 , const Block& M3) {
			typename Block::iterator p1;
			typename Block::const_iterator p3;
			for (p1 = M1.begin(),p3=M3.begin(); p3 != M3.end();p1++,p3++) {
				M2.apply(*p1,*p3);
			}
		}
	    
	
		/// User Left and Right vectors 
		void init (const Block& uu, const Block& vv) {
			casenumber = 1;
			u = uu;  //cout<<" u : ";u.write()<<endl;
			v = vv;  //cout<<" v : ";v.write()<<endl;
			_row = uu.rowdim();
			_col = vv.coldim();
			_value = Matrix(_field,u.rowdim(),v.coldim());
			Mul(_value,u,v);
		}

		// Random Left Matrix and Right Matrix
	
		void init (int n, int m) {
			casenumber = 1;
			_row = n;
			_col = m;
			typename Field::RandIter r(_field);

			u = Block(n , Vector(_BB->rowdim(),Zero));
			v = Block(m , Vector(_BB->coldim(),Zero));

			for (long i=0;i<n;i++)
				for (long j=0;j<_BB->rowdim();j++)
					r.random(u[i][j]);

			for (long i=0;i<m;i++)
				for (long j=0;j<_BB->coldim();j++)
					r.random(v[i][j]);

			_value = Matrix(_field,n,m);
			Mul(_value,u,v);
	  
		}

      

	};
 
}

#endif // __BLACKBOX_BLOCK_CONTAINER_BASE_H
