/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/switch/boolean.h
 * Copyright (C) 1999-2001 William J Turner
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *
 * ------------------------------------
 * 2002-08-20  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Brought this file into the current Linbox framework:
 *   - Renamed file as boolean.h
 *   - Renamed class boolean_switch as BooleanSwitch
 *   - Reindent
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BOOLEAN_H
#define __BOOLEAN_H

#include <vector>

namespace LinBox
{

/** Boolean switch object.
 * This is a switch predicate object that is applied
 * to two references to elements to switch them as needed 
 * by the \Ref{Butterfly Switching Network BlackBox Matrix Object}.
 */
class BooleanSwitch
{
    public:

	/** Constructor from an STL vector of booleans.
	 * The switch is applied using the vector of booleans.
	 * A true value means to swap the two elements, and a false
	 * value means not to.
	 * The apply function starts at the beginning of the vector moving 
	 * forward through it, and applyTranspose function starts at the end
	 * moving backwards.  Both repeat the vector after they pass through it.
	 * @param switches vector of switches
	 */
	BooleanSwitch (const std::vector<bool>& switches);

	/** Destructor.
	 */
	~BooleanSwitch (void) {}

	/** Apply switch function.
	 * Switches the elements in references according to current boolean
	 * value.  Swaps the elements if boolean is true, otherwise does nothing.
	 * It is templatized by the element type to be swapped.
	 * @return bool true if swapped, false otherwise
	 * @param x reference to first element to be switched
	 * @param y reference to second element to be switched
	 */
	template <class Element>
	bool apply (Element& x, Element& y) const;

	/** Apply switch transpose function.
	 * Switches the elements in references according to current boolean
	 * value.  Swaps the elements if boolean is true, otherwise does nothing.
	 * It is templatized by the element type to be swapped.
	 * @return bool true if swapped, false otherwise
	 * @param x reference to first element to be switched
	 * @param y reference to second element to be switched
	 */
	template <class Element>
	bool applyTranspose (Element& x, Element& y) const;

    private:

	// STL vector of boolean flags for switches
	std::vector<bool> _switches;

	// STL vector iterator and reverse iterator pointing to current switch
	// and its transpose
	mutable std::vector<bool>::const_iterator _iter;
	mutable std::vector<bool>::const_reverse_iterator _riter;

}; // class boolean_switch

inline BooleanSwitch::BooleanSwitch (const std::vector<bool>& switches)
	: _switches (switches)
{ 
	_iter = _switches.begin (); 
	_riter = _switches.rbegin (); 
}

template <class Element> 
inline bool BooleanSwitch::apply (Element& x, Element& y) const
{
	// If at end of vector, extend it
	if (_iter == _switches.end ()) {
		if (_switches.size () == 0)
			return false;
		else   
			_iter = _switches.begin ();
	}

	// If true value, swap
	if (*_iter)
		std::swap (x, y);

	// return with flag of *_iter and update _iter
	return *_iter++;
}
  
template <class Element> 
inline bool BooleanSwitch::applyTranspose (Element& x, Element& y) const
{
	// If at end of vector, extend it
	if (_riter == _switches.rend ()) {
		if (_switches.size () == 0)
			return false;
		else
			_riter = _switches.rbegin ();
	}

	// If true value, swap
	if (*_riter)
		std::swap (x, y);

	// return with flag of *_iter and update _iter
	return *_riter++;
}

}

#endif // __BOOLEAN_H
