/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/switch/cekstv.h
 * Copyright (C) 1999-2001 William J Turner
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *
 * ------------------------------------
 * 2002-08-20  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Brought this file into the current Linbox framework:
 *   - Renamed file as cekstv.h
 *   - Renamed class cekstv_switch as CekstvSwitch
 *   - Reindent
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __CEKSTV_H
#define __CEKSTV_H

#include <vector>

namespace LinBox
{

/** Butterfly switch object from preconditioner paper.
 * This is a switch predicate object that is applied
 * to two references to elements to switch them as needed 
 * by the \Ref{Butterfly Switching Network BlackBox Matrix Object}
 * following the exchange matrix introduced in "Efficient Matrix
 * Preconditioners for Black Box Linear Algebra" by Chen, Eberly, 
 * Kaltofen, Saunders, Turner, and Villard.
 * This class is templatized by the field in which the arithmetic
 * is done.
 */
template <class Field>
class CekstvSwitch
{
    public:

	/// Typedef
	typedef typename Field::Element Element;

	/** Constructor from a field and random field element generator.
	 * The switch is applied using the vector of field elements.
	 * If the current switch is marked by the field element a,
	 * and the two elements are x and y, the output from the switch
	 * is x' = x + a*y and y' = y + (x + a*y).
	 * The generator is used to create random field elements for setting the 
	 * switches.  Both apply function and applyTranspose function 
	 * use the random elements in the order they are generated.  Neither 
	 * can re-use random elements.  This means each time the matrix (or
	 * its transpose) is applied it will be a different matrix.
	 * @param F field in which arithmetic is done
	 * @param R random field element generator
	 */
	CekstvSwitch (const Field& F, const typename Field::RandIter& R);

	/** Constructor from a field and STL vector of field elements.
	 * The switch is applied using the vector of field elements.
	 * If the current switch is marked by the field element a,
	 * and the two elements are x and y, the output from the switch
	 * is x' = x + a*y and y' = y + (x + a*y).
	 * The apply function starts at the beginning of the vector moving 
	 * forward through it, and applyTranspose function starts at the end
	 * moving backwards.  Both repeat the vector after they pass through it.
	 * @param F field in which arithmetic is done
	 * @param switches vector of switches
	 */
	CekstvSwitch (const Field& F, const std::vector<Element>& switches);

	/** Copy constructor
	 */
	CekstvSwitch (const CekstvSwitch &s);

	/** Destructor.
	 */
	~CekstvSwitch () {}

	/** Apply switch function.
	 * Switches the elements in references according to the
	 * exchange matrix introduced in "Efficient Matrix
	 * Preconditioners for Black Box Linear Algebra" by Chen, Eberly, 
	 * Kaltofen, Saunders, Turner, and Villard and the current field element
	 * specified in the switch object.
	 * @return bool true if swapped, false otherwise
	 * @param x reference to first element to be switched
	 * @param y reference to second element to be switched
	 */
	bool apply (Element& x, Element& y) const;

	/** Apply switch transpose function.
	 * Switches the elements in references according to the
	 * transpose of the exchange matrix introduced in "Efficient Matrix
	 * Preconditioners for Black Box Linear Algebra" by Chen, Eberly, 
	 * Kaltofen, Saunders, Turner, and Villard and the current field element
	 * specified in the switch object.
	 * @return bool true if swapped, false otherwise
	 * @param x reference to first element to be switched
	 * @param y reference to second element to be switched
	 */
	bool applyTranspose (Element& x, Element& y) const;

   private:

	// Field in which arithemetic is done
	Field _F;

	// Random field element generator;
	mutable typename Field::RandIter _R;

	// STL vector of boolean flags for switches
	std::vector<Element> _switches;

	// STL vector iterator and reverse iterator pointing to current switch
	// and its transpose
	mutable typename std::vector<Element>::const_iterator _iter;
	mutable typename std::vector<Element>::const_reverse_iterator _riter;

	// temporary field element used in arithmetic
	mutable Element _temp;
};

template <class Field>
inline CekstvSwitch<Field>::CekstvSwitch (const Field& F, const typename Field::RandIter& R)
	: _F (F), _R (R)
{
	_iter = _switches.begin (); 
	_riter = _switches.rbegin (); 
}

template <class Field>
inline CekstvSwitch<Field>::CekstvSwitch (const Field& F, const std::vector<typename Field::Element>& switches)
	: _F (F), _R (_F), _switches (switches)
{ 
	_iter = _switches.begin (); 
	_riter = _switches.rbegin (); 
}

template <class Field>
inline CekstvSwitch<Field>::CekstvSwitch (const CekstvSwitch &s) 
	: _F (s._F), _R (s._R), _switches (s._switches)
{
	_iter = _switches.begin (); 
	_riter = _switches.rbegin (); 
}

template <class Field> 
inline bool CekstvSwitch<Field>::apply (typename Field::Element& x, typename Field::Element& y) const
{
	if (_switches.empty ()) {
		_R.random (_temp);
		_F.addin (x, _F.mulin (_temp, y));
		_F.addin (y, x);
	} else {    
		// If at end of vector, repeat it
		if (_iter == _switches.end ())
			_iter = _switches.begin ();

		_F.addin (x, _F.mul (_temp, *_iter++, y));
		_F.addin (y, x);
	}

	// return with swap flag
	return true;
}
  
template <class Field> 
inline bool CekstvSwitch<Field>::applyTranspose (typename Field::Element& x, typename Field::Element& y) const
{
	if (_switches.empty ()) {
		_F.addin (x, y);
		_R.random (_temp);
		_F.mulin (_temp, x);
		_F.addin (y, _temp);
	} else {    
		// If at end of vector, extend it
		if (_riter == _switches.rend ())
			_riter = _switches.rbegin ();

		_F.addin (x, y);
		_F.addin (y, _F.mul (_temp, *_riter++, x));
	}

	// return with swap flag
	return true;
}

}

#endif // __CEKSTV_H
