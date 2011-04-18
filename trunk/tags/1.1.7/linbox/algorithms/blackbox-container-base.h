/* linbox/algorithms/blackbox-container.h
 * Copyright (C) 1999, 2001 Jean-Guillaume Dumas
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
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

#ifndef __LINBOX_blackbox_container_base_H
#define __LINBOX_blackbox_container_base_H

// ================================================================
// Base ForwardIterator wrapper for BlackBoxes
// Have to be provided :
// - launch : launches the following computation
// - wait   : waits for the end of the current computation
// ================================================================


#include "linbox/vector/vector-domain.h"

namespace LinBox 
{

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif

/** \brief A base class for BlackboxContainer.
  * The primary member function is begin().

  * It returns an iterator which after i increments (++) dereferences to 
  * $v^T A^i u$, for $v$ and $u$ determined by the form of construction.
  * It is designed to be used with implementations of Berlekamp-Massey
  * such as MasseyDom.
  *
  * Subclasses complete the implementation by defining _launch() and _wait().
  */

template<class Field, class Blackbox>
class BlackboxContainerBase {
    public:
	typedef typename Field::Element Element;      

        //-- Constructors
	BlackboxContainerBase () {} 

	BlackboxContainerBase (const Blackbox *BB, const Field &F)
		: _F (F), _VD (F), _BB (BB), _size (MIN (BB->rowdim (), BB->coldim ())) { _size <<= 1; }
	
	// Pascal Giorgi 16.02.2004
	BlackboxContainerBase (const Blackbox *BB, const Field &F, unsigned long size)
		: _F (F), _VD (F), _BB (BB), _size (size) {}

	virtual ~BlackboxContainerBase ()
	{ }//delete _BB; }

	class const_iterator {
		BlackboxContainerBase<Field, Blackbox> &_c;
	public:
		//const_iterator () {} // BB ??
		const_iterator (BlackboxContainerBase<Field, Blackbox> &C) : _c (C) {}
		const_iterator &operator ++ () { _c._launch (); return *this; }
		const Element  &operator *  () { _c._wait ();   return _c.getvalue (); }
	};

	const_iterator begin () { return const_iterator (*this); }
	const_iterator end   () { return const_iterator (); }

	long         size     () const { return _size; }
	const Field &getField () const { return _F; }
	Blackbox    *getBB    () const { return _BB; }

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

	Field                _F;
	VectorDomain<Field>  _VD;
	const Blackbox            *_BB;
    
	long                 _size;

        // BDS 22.03.03
	long                 casenumber;
	std::vector<Element>    u, v;
	Element              _value;

	const Element &getvalue() { return _value; }

	//-------------- 
	/// Initializers
	//--------------  

        /// User Left and Right vectors 
	template<class Vector1, class Vector2>
	Element &init (const Vector1& uu, const Vector2& vv) {
		casenumber = 1;
		u.resize(uu.size());
		std::copy(uu.begin(),uu.end(),u.begin());
		//u = uu;
		v.resize(vv.size());
		std::copy(vv.begin(),vv.end(),v.begin());
		//v = vv;
                    // JGD 22.03.03
// 		return _VD.dot (_value, u, u);
		return _VD.dot (_value, u, v);
	}

        /// Random Left vectors, Zero Right vector
	template<class RandIter>
	Element &init (RandIter& g)
	{
		casenumber = 1;
		u.resize (_BB->coldim ());
		for (long i = u.size (); i--;)
			g.random (u[i]);
		v.resize (_BB->rowdim ());
		return _VD.dot (_value, u, u);
	}

        /// User Left vectors, Zero Right vector
	template<class Vector>
	Element &init (const Vector& uu) {
		casenumber = 1;
		u.resize(uu.size());
		std::copy(uu.begin,uu.end(),u.begin());
		v.resize (_BB->rowdim ());
		return _VD.dot (_value, u, u);
	}
};
 
}

#endif // __LINBOX_blackbox_container_base_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
