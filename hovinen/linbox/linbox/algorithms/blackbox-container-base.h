/* -*- mode: c; style: linux -*- */

/* linbox/src/library/objects/algorithms/blackbox/blackbox-container.h
 * Copyright (C) 1999, 2001 Jean-Guillaume Dumas, Bradford Hovinen
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
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

// ================================================================
// LinBox Project 1999
// Base ForwardIterator wrapper for BlackBoxes
// Have to be provided :
// - launch : launches the following computation
// - wait   : waits for the end of the current computation
// Time-stamp: <26 May 00 17:38:58 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================

#ifndef __BLACKBOX_CONTAINER_BASE_H
#define __BLACKBOX_CONTAINER_BASE_H

#include "linbox/blackbox/archetype.h"

namespace LinBox 
{

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

/** BlackboxContainerBase is a base class for BlackboxContainer
  * begin() is the primary member function.
  * It returns an iterator which after i increments (++) dereferences to 
  * $v^T A^i u$, for $v$ and $u$ determined by the form of construction.
  * It is designed to be used with implementations of Berlekamp-Massey
  * such as MasseyDom.
  *
  * Subclasses complete the implementation by defining _launch() and _wait().
  */

template<class Field, class Vector>
class BlackboxContainerBase {
    public:
	typedef Blackbox_archetype<Vector> Blackbox;
	typedef typename Field::Element Element;

        //-- Constructors
	BlackboxContainerBase () {} 

	BlackboxContainerBase (const Blackbox *BD, const Field &F)
		: _field (F), _BB (BD->clone ()), _size (MIN (BD->rowdim (), BD->coldim ()) << 1) {}
    
	class const_iterator {
		BlackboxContainerBase<Field, Vector> &_c;
	public:
		const_iterator () {}
		const_iterator (BlackboxContainerBase<Field, Vector> &C) : _c (C) {}
		const_iterator &operator ++ () { _c._launch (); return *this; }
		const Element &operator * () { _c._wait (); return _c.getvalue(); }
	};

	const_iterator begin () { return const_iterator (*this); }
	const_iterator end () { return const_iterator (); }

	long size () { return _size; }
	const Field &getField () const { return _field; }
	Blackbox *getBB () const { return _BB; }

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
	Blackbox *_BB;
    
	long _size;

	bool even;
	Vector u, v;
	Element _value;
	const Element &getvalue() { return _value; }

//-------------- 
/// Initializers
//--------------  
        /// User Left and Right vectors 
	Element &init (const Vector& uu, const Vector& vv) {
		even = 1;
		u = uu;
		v = vv;
		return DOTPROD (_value, u, u);
	}

        /// Random Left vectors, Zero Right vector
	template<class RandIter>
	Element &init (RandIter& g) {
		even = 1;
		u.resize (_BB->coldim ());
		for (long i = u.size (); i--;)
			g.random (u[i]);
		v.resize (_BB->rowdim ());
		return DOTPROD (_value, u, u);
	}

        /// User Left vectors, Zero Right vector
	Element &init (const Vector& uu) {
		even = 1;
		u = uu;
		v.resize (_BB->rowdim ());
		return DOTPROD (_value, u, u);
	}

//-------------- 
/// Operators
//--------------  
        /// Generic dot product using the container domain
	template <class A1, class A2>
	Element &DOTPROD (Element &coeff, const A1 &u, const A2 &v) {
		_field.mul (coeff, u[0], v[0]);
		for (long k = u.size () - 1; k > 0; --k)
			_field.axpyin (coeff, u[k], v[k]);
		return coeff;
	}

        /// Generic axpy using the container domain
	template <class A1, class A2, class A3>
	A1 &AXPY (A1 &u, const Element &coeff, const A2 &v, const A3 &w) {
		for (long k = u.size() - 1; k >= 0; --k)
			_field.axpy (u[k], coeff, v[k], w[k]);
		return u;
	}

        /// u <-- u + c * v
	template <class A1, class A2>
	A1 &AXPYIN (A1 &u, const Element &coeff, const A2 &v) {
		for (long k = u.size () - 1 ; k >= 0; --k)
			_field.axpyin (u[k], coeff, v[k]);
		return u;
	}

        /// u <-- u - c * v
	template <class A1, class A2>
	A1& AXMYIN (A1& u, const Element &coeff, const A2& v) {
		for (long k = u.size () - 1 ; k >= 0; --k)
			_field.axmyin (u[k], coeff, v[k]);
		return u;
	}

        /// u <-- u + v
	template <class A1, class A2>
	A1& ADDIN (A1& u, const A2& v) {
		for (long k = u.size () - 1 ; k >= 0; --k)
			_field.addin (u[k], v[k]);
		return u;
	}

        /// u <-- c * u
	template <class A1>
	A1 &MULIN (A1& u, const Element& coeff) {
		for (long k = u.size () - 1 ; k >= 0; --k)
			_field.mulin (u[k], coeff);
		return u;
	}

        /// u <-- c * v
	template <class A1, class A2>
	A1 &MUL (A1& u, const Element &coeff, const A2 &v) {
		for (long k = u.size () - 1 ; k >= 0; --k)
			_field.mul (u[k], coeff, v[k]);
		return u;
	}

        /// u <-- c * u + v
	template <class A1, class A2>
	A1& AXINPY (A1 &u, const Element &coeff, const A2 &v) {
		Element tmp;
		for (long k = u.size () - 1 ; k >= 0; --k) {
			_field.axpy (tmp, u[k], coeff, v[k]);
			_field.assign (u[k], tmp);
		}
		return u;
	}

        /// u <-- c * u - v
	template <class A1, class A2>
	A1 &AXINMY (A1 &u, const Element &coeff, const A2 &v) {
		Element tmp;
		for (long k = u.size () - 1 ; k >= 0; --k) {
			_field.axmy (tmp, u[k], coeff, v[k]);
			_field.assign (u[k], tmp);
		}
		return u;
	}
};
 
}

#endif // __BLACKBOX_CONTAINER_BASE_H
